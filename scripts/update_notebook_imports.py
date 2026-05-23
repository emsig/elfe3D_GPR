import json
from pathlib import Path

root = Path('examples')
notebooks = list(root.rglob('*.ipynb'))
changed = []
for nb_path in notebooks:
    data = json.loads(nb_path.read_text(encoding='utf-8'))
    modified = False
    for cell in data.get('cells', []):
        if cell.get('cell_type') != 'code':
            continue
        source = cell.get('source', [])
        if isinstance(source, str):
            lines = source.splitlines(keepends=True)
        else:
            lines = source
        new_lines = []
        for line in lines:
            if "sys.path.insert(0, os.path.join(MASTER, 'io'))" in line or \
               "sys.path.insert(0, os.path.join(MASTER, 'io', 'inputs'))" in line or \
               "sys.path.insert(0, os.path.join(MASTER, 'io', 'outputs'))" in line:
                modified = True
                continue
            if line.strip().startswith('from survey import '):
                modified = True
                new_lines.append(line.replace('from survey import ', 'from elfe3d_gpr.inputs.survey import '))
                continue
            if line.strip().startswith('from runner import '):
                modified = True
                new_lines.append(line.replace('from runner import ', 'from elfe3d_gpr.runner import '))
                continue
            if line.strip().startswith('from anomalies import '):
                modified = True
                new_lines.append(line.replace('from anomalies import ', 'from elfe3d_gpr.inputs.anomalies import '))
                continue
            if line.strip().startswith('from fieldreader import '):
                modified = True
                new_lines.append(line.replace('from fieldreader import ', 'from elfe3d_gpr.outputs.fieldreader import '))
                continue
            if line.strip().startswith('from postprocess import '):
                modified = True
                new_lines.append(line.replace('from postprocess import ', 'from elfe3d_gpr.outputs.postprocess import '))
                continue
            if line.strip().startswith('from visualize import '):
                modified = True
                new_lines.append(line.replace('from visualize import ', 'from elfe3d_gpr.outputs.visualize import '))
                continue
            if line.strip().startswith('from materials import '):
                modified = True
                new_lines.append(line.replace('from materials import ', 'from elfe3d_gpr.inputs.materials import '))
                continue
            if line.strip().startswith('from receivers import '):
                modified = True
                new_lines.append(line.replace('from receivers import ', 'from elfe3d_gpr.inputs.receivers import '))
                continue
            if line.strip().startswith('from sources import '):
                modified = True
                new_lines.append(line.replace('from sources import ', 'from elfe3d_gpr.inputs.sources import '))
                continue
            if line.strip().startswith('from solver import '):
                modified = True
                new_lines.append(line.replace('from solver import ', 'from elfe3d_gpr.inputs.solver import '))
                continue
            if line.strip().startswith('from pml import '):
                modified = True
                new_lines.append(line.replace('from pml import ', 'from elfe3d_gpr.inputs.pml import '))
                continue
            if line.strip().startswith('from writetetgenpoly import '):
                modified = True
                new_lines.append(line.replace('from writetetgenpoly import ', 'from elfe3d_gpr.inputs.writetetgenpoly import '))
                continue
            if line.strip().startswith('from writeinputfiles import '):
                modified = True
                new_lines.append(line.replace('from writeinputfiles import ', 'from elfe3d_gpr.inputs.writeinputfiles import '))
                continue
            new_lines.append(line)
        if modified:
            cell['source'] = ''.join(new_lines) if isinstance(source, str) else new_lines
    if modified:
        nb_path.write_text(json.dumps(data, indent=1, ensure_ascii=False), encoding='utf-8')
        changed.append(str(nb_path))
print('modified notebooks:')
for path in changed:
    print(path)
