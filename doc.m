%% geoid
%
% Function to generate a n-level geoid shape, defined by its vertices
% set and the corresponding triangulation. Vertices belong to the unit sphere.
%
% Algorithm principle is based on sampling the 20 icosahedron triangles.
%
% Author & support : nicolas.douillet (at) free.fr, 2017-2021.
%
%% Syntax
% [V, T] = geoid(n)
%
% [V, T] = geoid(n, option_display)
%
%% Description
% [V, T] = geoid(n) computes n-level geoid vertices coordinates V
% and its corresponding triangulation T.
%
% [V, T] = geoid(n, option_display) display the resulting geoid if
% the boolean option_display is either logical or numeric true.
%
%% See also 
%
% <https://fr.mathworks.com/matlabcentral/fileexchange/64284-sphere-homeomorphic-surface-quasi-isotropic-sampling?s_tid=prof_contriblnk sphere_homeo_sfc_isotropic_splg> |
% <https://fr.mathworks.com/help/matlab/ref/sphere.html?s_tid=srchtitle sphere> |
% <https://fr.mathworks.com/matlabcentral/fileexchange/64395-sample_triangle?s_tid=prof_contriblnk sample_triangle>
%
%% Input arguments
%
% - n : positive integer scalar double (n > 0), the level which corresponds to the
%       number of steps to sample the initial triangles in.
%
% - option_display : logical *true (1) / false (0).
%
%% Output arguments
%
%        [ |  |  |]
% - V = [Vx Vy Vz], real matrix double, the geoid vertices X, Y, Z coordinates. Size(V) = [nb_vertices,3], with nb_vertices < (20*(n+1)*(n+2)/2.
%        [ |  |  |]       
%
%        [ |  |  |]
% - T = [i0 i1 i2], positive integer matrix double, the triangulation. Size(T) = [nb_triangles,3], with nb_triangles < 20*n^2.
%        [ |  |  |]
%
%% Example #1
%
% Basic icosahedron (level = 1)
n = 1;
option_display = true;
[V, T] = geoid(n, option_display);

%% Example #2
%
% 6th level geoid
n = 6;
option_display = true;
[V, T] = geoid(n, option_display);