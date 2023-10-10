
function [ x0, y0, z0, th_rfl ] = make_reflec( spl2d, vec_a, vec_v, ui, vj )

% function [ x0, y0, z0, th_rfl ] = make_reflec( spl2d, vec_a, vec_v, ui, vj  )
%
% compute reflection angles theta [deg] for light source direction vec_a and view direction vec_v
% sampling the spline at u-parameters ui and v-parameters vj
%
% vec_a, vec_v can be [x,y,z] or [az,el], with az,el in degrees

rad = pi / 180;

if (length(vec_a)==2)
  a_az = vec_a(1);
  a_el = vec_a(2);
  vec_a = [ cos(a_az*rad)*cos(a_el*rad), sin(a_az*rad)*cos(a_el*rad), sin(a_el*rad) ]';
end
if (length(vec_v)==2)
  v_az = vec_v(1);
  v_el = vec_v(2);
  vec_v = [ cos(v_az*rad)*cos(v_el*rad), sin(v_az*rad)*cos(v_el*rad), sin(v_el*rad) ]';
end
if (size(vec_a,2)>1)
   vec_a = vec_a';
end
if (size(vec_v,2)>1)
   vec_v = vec_v';
end

% evaluate spline surface at positions ui x vj

idebug = 0;
[~, x0, y0, z0] = eval_2dspline(spl2d, ui, vj, [], idebug);

% construct first normalized tangent vector at each ui x vj

[dx1, dy1, dz1] = eval_2dspline_deriv(spl2d, ui, vj, 1);
nrm = sqrt(dx1.^2 + dy1.^2 + dz1.^2);
t1(:,:,1) = dx1./nrm; t1(:,:,2) = dy1./nrm; t1(:,:,3) = dz1./nrm;

% construct second normalized tangent vector at each ui x vj

[dx2, dy2, dz2] = eval_2dspline_deriv(spl2d, ui, vj, 2);
nrm = sqrt(dx2.^2 + dy2.^2 + dz2.^2);
t2(:,:,1) = dx2./nrm; t2(:,:,2) = dy2./nrm; t2(:,:,3) = dz2./nrm;

% construct normalized normal vector at each ui x vj

nvec = cross(t1, t2);
nrm = sqrt(sum(nvec.^2, 3));
nvec(:,:,1) = nvec(:,:,1) ./ nrm;
nvec(:,:,2) = nvec(:,:,2) ./ nrm;
nvec(:,:,3) = nvec(:,:,3) ./ nrm;

% compute projection of v onto plane P

v_proj = vec_v - (vec_a'*vec_v)*vec_a; v_proj = v_proj / norm(v_proj);
v_orth = cross(vec_a, v_proj);

% compute theta at surface evaluation positions

nu = length(ui);
nv = length(vj);
th_rfl = zeros(nu, nv);

if (0==1)

   for iu = 1 : nu
      for jv = 1 : nv
         n_i   = squeeze(nvec(iu,jv,:));
         v_t   = vec_v - (n_i'*vec_v)*n_i;
         vec_r = 2 * v_t - vec_v;

         r_proj = vec_r - (vec_a'*vec_r)*vec_a; r_proj = r_proj / norm(r_proj);
         r_p   = (-r_proj'*v_proj);
         r_o   = (-r_proj'*v_orth);
         th_rfl(iu,jv) = atan2(r_o, r_p) / rad;
      end
   end

else

   one  = ones(nu,nv);

   % n_i   = squeeze(nvec(iu,jv,:));
   % v_t   = vec_v - (n_i'*vec_v) * n_i;
   % vec_r = 2 * v_t - vec_v;

   dot_nv = zeros(nu,nv);
   for k = 1 : 3
      dot_nv = dot_nv + vec_v(k) * nvec(:,:,k);
   end

   vec_r  = zeros(nu,nv,3);
   dot_ar = zeros(nu,nv);
   for k = 1 : 3
      vec_r(:,:,k) = vec_v(k) - 2 * dot_nv .* nvec(:,:,k);
      dot_ar = dot_ar + vec_a(k) * vec_r(:,:,k);
   end

   % r_proj = vec_r - (vec_a'*vec_r)*vec_a;
   % r_p   = (-r_proj'*v_proj);
   % r_o   = (-r_proj'*v_orth);

   r_proj = zeros(nu,nv,3);
   for k = 1 : 3
      r_proj(:,:,k) = vec_r(:,:,k) - dot_ar * vec_a(k);
   end

   dot_rvp = zeros(nu,nv);
   dot_rvo = zeros(nu,nv);
   for k = 1 : 3
      dot_rvp = dot_rvp - r_proj(:,:,k) * v_proj(k);
      dot_rvo = dot_rvp - r_proj(:,:,k) * v_orth(k);
   end

   th_rfl = atan2(dot_rvo, dot_rvp) / rad;

end

end % make_reflec
