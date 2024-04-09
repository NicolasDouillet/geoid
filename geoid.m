function [V, T] = geoid(n, option_display)
% geoid : function to generate a n-level geoid shape, defined by its
% vertices set and the corresponding triangulation. Vertices belong to the unit sphere.
% Algorithm principle is based on sampling the 20 triangles of the icosahedron.
%
% Author : nicolas.douillet (at) free.fr, 2017-2024.
%
%
% Syntax
% [V, T] = geoid(n)
% [V, T] = geoid(n, option_display)
%
%
% Description
% [V, T] = geoid(n) computes n-level geoid vertices coordinates V
% and its corresponding triangulation T.
%
% [V, T] = geoid(n, option_display) display the resulting geoid if
% the boolean option_display is is either logical or numeric true.
%
%
% See also : SPHERE
%
%
% Input arguments
%
% - n : positive integer scalar double (n > 0), the level which corresponds to the
%       number of steps to sample the initial triangles in.
%
% - option_display : logical *true (1) / false (0).
%
%
% Output arguments
%
%       [ |  |  |]
% - V = [Vx Vy Vz], real matrix double, the geoid vertices X, Y, Z coordinates. Size(V) = [nb_vertices,3], with nb_vertices < (20*(n+1)*(n+2)/2.
%       [ |  |  |]       
%
%       [ |  |  |]
% - T = [i0 i1 i2], positive integer matrix double, the triangulation. Size(T) = [nb_triangles,3], with nb_triangles < 20*n^2.
%       [ |  |  |]
%
%
% Example #1
%
% Basic icosahedron (level = 1)
% n = 1;
% option_display = true;
% [V, T] = geoid(n, option_display);
%
%
% Example #2
%
% 6th level geoid
% n = 6;
% option_display = true;
% [V, T] = geoid(n, option_display);


% Input parsing
if nargin < 2
    option_display = true;
end

assert(nargin >= 1,'Error : not enough input arguments.');
assert(nargin < 3,'Error : too many input arguments.');
assert(isreal(n) && n >= 1 && rem(n,1) == 0,'Error : level must be a stricly positive integer.');


% Body
% Step I : build icosahedron
Rho = 1;
phi_n = 0.5*(1+sqrt(5));
centre_angle = 2*asin(1/Rho/sqrt(phi_n*sqrt(5)));

Rmz = @(theta)[cos(theta) -sin(theta) 0;...
               sin(theta)  cos(theta) 0;...
               0          0           1];

% 1st equilateral triangle
V1 = [0 0 Rho]';
V2 = [Rho*sin(centre_angle) 0 Rho*cos(centre_angle)]';
V3 = Rmz(0.4*pi)*V2;

% Lower base triangle with /O symetry
V12 = -V1;
V11 = -V2;
V10 = -V3;

% (12) vertices set coordinates vector
V4 = Rmz(0.4*pi)*V3;
V5 = Rmz(0.8*pi)*V3;
V6 = Rmz(1.2*pi)*V3;
V9 = Rmz(0.4*pi)*V10;
V8 = Rmz(0.8*pi)*V10;
V7 = Rmz(1.2*pi)*V10;

V = [V1 V2 V3 V4 V5 V6 V7 V8 V9 V10 V11 V12]';

T = [1 2 3; % top pentagone (5 triangles)
     1 3 4;
     1 4 5;
     1 5 6;
     1 6 2; ...
     
     12 11 10; % bottom pentagone (5 triangles)
     12 10 9;
     12 9 8;
     12 8 7;
     12 7 11; ...
     
     2 3 8; % centre belt (10 triangles)
     3 8 7;
     3 4 7;
     4 7 11;
     4 5 11;
     5 11 10;
     5 6 10;
     6 10 9;
     6 2 9;
     2 9 8];

 
% Step II : oversample triangles by creating new vertices ; link vertices to create new triangles
nt = n^2; % nb new triangles
nv = (n+1)*(n+2)/2; % nb new vertices

T_new = zeros(20*nt,3);
V_new = zeros(10*(n+1)*(n+2),3);

for j = 1:size(T,1)
    
    [new_sub_V, new_sub_T] = sample_triangle(V(T(j,1),:)', V(T(j,2),:)', V(T(j,3),:)', n);
    new_sub_T = new_sub_T + (j-1)*nv; % update triangle indices
    
    T_new((j-1)*nt+1:j*nt,:) = new_sub_T;
    V_new((j-1)*nv+1:j*nv,:) = new_sub_V;
    
end

T = T_new;


% Step III : remove duplicated vertices and re-index
V_unique = zeros(size(unique(V_new,'rows'),1),3);
curr_vtx_idx = 0;
skiplist = [];

for k = 1:size(V_new,1)
    
    if isempty(find(skiplist == k, 1))
        
        curr_vtx_idx = curr_vtx_idx + 1;
        V_unique(curr_vtx_idx,:) = V_new(k,:);
        f = find(all(bsxfun(@eq,V_new,V_new(k,:)),2));
        
        if size(f,1) > 1 % at least one duplication
            
            skiplist = [skiplist; f(2:end,1)];
            
        end
        
        for s = 1:size(f,1)
            
            g = (T == f(s,1));
            T(g) = curr_vtx_idx; % Replace corresponding indices in T triangles list
            
        end
        
    end
    
end

V = V_unique;


% Step IV : vertices normalization / projection on the sphere surface
norm_V = sqrt(sum(V.^2,2));
V = V./repmat(norm_V,[1,3]);


% Remove duplicated vertices
[V, T] = remove_duplicated_vertices(V, T);


% Remove duplicated triangles
T = unique(sort(T,2), 'rows', 'stable');


% Step V (optional) : display
if option_display
    
    figure;
    set(gcf, 'Color', [0 0 0]);
    TRI = triangulation(T, V(:,1), V(:,2), V(:,3));
    t = trimesh(TRI, 'FaceAlpha',0.2);
    colormap([0 1 0]);
    t.FaceColor = [0 1 0];
    t.LineWidth = 2;
    set(gca, 'Color', [0 0 0], 'XColor', [1 1 1], 'YColor', [1 1 1], 'ZColor', [1 1 1]);
    axis equal;
    
end

end % geoid


% sample_triangle subfunction
function [V, T] = sample_triangle(V1, V2, V3, nbstep)


% Create sampling grid
global Ndim;

Ndim = size(V1,1);

% (V1V2, V1V3) base
u = (V2 - V1);
v = (V3 - V1);

V = zeros(sum(1:nbstep+1),Ndim);

nu = u / norm(u);
nv = v / norm(v);
stepu = norm(u) / nbstep;
stepv = norm(v) / nbstep;
k = 1;

% Sampling & vertices generation
for m = 0:nbstep
    
    for n = 0:nbstep
        
        if m+n <= nbstep % in (V1,V2,V3) triangle conditions ; indices # nb segments
            
            % translation vector
            tv = m*stepu*nu + n*stepv*nv;
            V(k,:) = (V1 + tv)';
            k = k+1;
            
        end
        
    end
    
end


% Index triplets list construction
T = zeros(nbstep^2,3);
row_length = 1 + nbstep;
cum_row_length = row_length;
row_idx = 1;
p = 1;

while p <= nbstep^2 && row_length > 1
    
     i = p;
    
    if p < 2 % "right" triangle serie only
        
        while (i < cum_row_length)
            
            T(row_idx,:) = [i i+1 i+row_length];
            row_idx = row_idx + 1;
            i = i +1;
            
        end
        
        row_length = row_length - 1;
        cum_row_length = cum_row_length + row_length;
        p = p + row_length+1;
        
    else
        
        % Since p >= 2
        while i < cum_row_length % both triangle series
            
            T(row_idx,:) = [i i+1 i+row_length];
            row_idx = row_idx + 1;            
            T(row_idx,:) = [i i-row_length i+1]; % + upside-down triangles serie
            row_idx = row_idx + 1;            
            i = i +1;
            
        end
        
        row_length = row_length - 1;
        cum_row_length = cum_row_length + row_length;
        p = p + row_length+1;
        
    end
    
end

T = sort(T,2);
T = unique(T,'rows','stable');


end % sample_triangle


% Remove duplicated vertices subfunction
function [V_out, T_out] = remove_duplicated_vertices(V_in, T_in)


tol = 1e4*eps;
[V_out,~,n] = uniquetol(V_in,tol,'ByRows',true);
T_out = n(T_in);


end % remove_duplicated_vertices