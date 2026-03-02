clear variables
clc
% plan view just below source level
z=0.01;
% compute the fields in the horizontal plane in a custom grid
xmin = 0.2;
xmax = 4.8;
nx = 256;
x = linspace(xmin, xmax, nx);
dx = x(2) - x(1);  % spacing
xx_pre = [0.0001, 0.001, 0.01];
xx_post = x(end) + dx : dx : 5.0;
xx = [xx_pre, x, xx_post];
xx = sort(unique([-fliplr(xx), xx]));
yy=xx;
[x,y]=meshgrid(xx,yy);
% horizontal distance
r=sqrt(x.^2+y.^2);
% 3D distance
R=sqrt(r.^2+z^2);
% free space parameters
c0=299792458;
mu=4e-7*pi;
% frequency
freq=1e8;
% parameters and fields
smu=2i*pi*freq*mu;
eta0=2i*pi*freq/(c0^2*mu);
grg0=sqrt(smu.*eta0);
exx=((3*(x./R).^2-1).*(1 + grg0*R) ...
         + ((x./R).^2-1).*(grg0*R).^2).*exp(-grg0*R)./(4*pi*eta0*R.^3);
eyx=x.*y.*(3 + 3*grg0*R + (grg0*R).^2).*exp(-grg0*R)./(4*pi*eta0*R.^5);
ezx=x.*z.*(3 + 3*grg0*R + (grg0*R).^2).*exp(-grg0*R)./(4*pi*eta0*R.^5);

% amplitude and phase plots of all three vector components
components = {exx, eyx, ezx};
titles = {'E_{xx}(x,y,z = ', 'E_{yx}(x,y,z = ', 'E_{zx}(x,y,z = '};
filenames = {['exx0-xy-', num2str(freq*1e-6,3), 'MHz-xmax', num2str(xmax), '.png'], ...
             ['eyx0-xy-', num2str(freq*1e-6,3), 'MHz-xmax', num2str(xmax), '.png'], ...
             ['ezx0-xy-', num2str(freq*1e-6,3), 'MHz-xmax', num2str(xmax), '.png']};

for i = 1:3
    % Create a figure with subplots for amplitude and phase
    figure('Position', [100, 100, 1200, 600]) % Set figure size for better resolution
    
    % Amplitude subplot
    subplot(1, 2, 1)
    imagesc(xx, yy, log10(abs(components{i})));
    axis([-xmax xmax -xmax xmax])
    xlabel('x-distance (m)')
    ylabel('y-distance (m)')
    title(['log_{10}|', titles{i}, num2str(z,3), ' m)| at f = ', num2str(freq*1e-6,3), ' MHz'], 'Interpreter', 'tex')
    set(gca,'fontsize',18,'plotboxaspectratio',[1 1 1])
    colorbar; % Add colorbar
    
    % Phase subplot
    subplot(1, 2, 2)
    imagesc(xx, yy, angle(components{i}));
    axis([-xmax xmax -xmax xmax])
    xlabel('x-distance (m)')
    ylabel('y-distance (m)')
    title(['\phi_{', titles{i}, num2str(z,3), ' m)} at f = ', num2str(freq*1e-6,3), ' MHz'], 'Interpreter', 'tex')
    set(gca,'fontsize',18,'plotboxaspectratio',[1 1 1])
    colorbar; % Add colorbar
    
    % Save the figure
    print('-dpng', '-r300', filenames{i}) % Save as PNG with high resolution
end

% write a csv file with all three components
csv_filename = ['whole_space_data_z', num2str(z, 3), '_', num2str(freq*1e-6, 3), 'MHz_xmax', num2str(xmax), '.csv'];
header = {'x', 'y', ...
    'E_xx_amplitude', 'E_xx_phase', 'E_xx_real', 'E_xx_imag', ...
    'E_yx_amplitude', 'E_yx_phase', 'E_yx_real', 'E_yx_imag', ...
    'E_zx_amplitude', 'E_zx_phase', 'E_zx_real', 'E_zx_imag'};
data = [x(:), y(:), ...
    abs(exx(:)), angle(exx(:)), real(exx(:)), imag(exx(:)), ...
    abs(eyx(:)), angle(eyx(:)), real(eyx(:)), imag(eyx(:)), ...
    abs(ezx(:)), angle(ezx(:)), real(ezx(:)), imag(ezx(:))];
% Open file and write header
output_data = [header; num2cell(data)];
writecell(output_data, csv_filename);

% plan view approximately three wavelengths below source level
z=1;
R=sqrt(r.^2+z^2);
exx=((3*(x./R).^2-1).*(1 + grg0*R) ...
         + ((x./R).^2-1).*(grg0*R).^2).*exp(-grg0*R)./(4*pi*eta0*R.^3);
eyx=x.*y.*(3 + 3*grg0*R + (grg0*R).^2).*exp(-grg0*R)./(4*pi*eta0*R.^5);
ezx=x.*z.*(3 + 3*grg0*R + (grg0*R).^2).*exp(-grg0*R)./(4*pi*eta0*R.^5);

% amplitude and phase plots of all three vector components
components = {exx, eyx, ezx};
titles = {'E_{xx}(x,y,z = ', 'E_{yx}(x,y,z = ', 'E_{zx}(x,y,z = '};
filenames = {['exx1-xy-', num2str(freq*1e-6,3), 'MHz-xmax', num2str(xmax), '.png'], ...
             ['eyx1-xy-', num2str(freq*1e-6,3), 'MHz-xmax', num2str(xmax), '.png'], ...
             ['ezx1-xy-', num2str(freq*1e-6,3), 'MHz-xmax', num2str(xmax), '.png']};

for i = 1:3
    % Create a figure with subplots for amplitude and phase
    figure('Position', [100, 100, 1200, 600]) % Set figure size for better resolution
    
    % Amplitude subplot
    subplot(1, 2, 1)
    imagesc(xx, yy, log10(abs(components{i})));
    axis([-xmax xmax -xmax xmax])
    xlabel('x-distance (m)')
    ylabel('y-distance (m)')
    title(['log_{10}|', titles{i}, num2str(z,3), ' m)| at f = ', num2str(freq*1e-6,3), ' MHz'], 'Interpreter', 'tex')
    set(gca,'fontsize',18,'plotboxaspectratio',[1 1 1])
    colorbar; % Add colorbar
    
    % Phase subplot
    subplot(1, 2, 2)
    imagesc(xx, yy, angle(components{i}));
    axis([-xmax xmax -xmax xmax])
    xlabel('x-distance (m)')
    ylabel('y-distance (m)')
    title(['\phi_{', titles{i}, num2str(z,3), ' m)} at f = ', num2str(freq*1e-6,3), ' MHz'], 'Interpreter', 'tex')
    set(gca,'fontsize',18,'plotboxaspectratio',[1 1 1])
    colorbar; % Add colorbar
    
    % Save the figure
    print('-dpng', '-r300', filenames{i}) % Save as PNG with high resolution
end
csv_filename = ['whole_space_data_z', num2str(z, 3), '_', num2str(freq*1e-6, 3), 'MHz_xmax', num2str(xmax), '.csv'];
header = {'x', 'y', ...
    'E_xx_amplitude', 'E_xx_phase', 'E_xx_real', 'E_xx_imag', ...
    'E_yx_amplitude', 'E_yx_phase', 'E_yx_real', 'E_yx_imag', ...
    'E_zx_amplitude', 'E_zx_phase', 'E_zx_real', 'E_zx_imag'};
data = [x(:), y(:), ...
    abs(exx(:)), angle(exx(:)), real(exx(:)), imag(exx(:)), ...
    abs(eyx(:)), angle(eyx(:)), real(eyx(:)), imag(eyx(:)), ...
    abs(ezx(:)), angle(ezx(:)), real(ezx(:)), imag(ezx(:))];

% Write header and data to a CSV file
output_data = [header; num2cell(data)];
writecell(output_data, csv_filename);
