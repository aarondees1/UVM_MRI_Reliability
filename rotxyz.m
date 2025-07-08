function a = rotxyz(b)
%ROTXYZ Computes the xyz-convention rotation transformation matrix.
%       ROTXYZ(B) returns the rotation matrix based on three Euler
%       angles in B.  B must be of length 3.  This convention is
%       also know as Tait-Bryan angles or the 321 sequence.
%
%       Reference:
%            Goldstein H: Classical Mechanics. 2nd Ed. Addison-Wesley
%            Publishing Co., Reading, Mass., 1980, pp 143-8, 608-10
%
%       31-Aug-94
%

% Check if input vector b has exactly 3 elements.
if (prod(size(b))==3)
  % Compute Z-axis rotation matrix (rotation by b(3)).
  c = cos(b(3)); % Cosine of Z-axis angle
  s = sin(b(3)); % Sine of Z-axis angle
  D = [c s 0; -s c 0; 0 0 1]; % Z-axis rotation matrix

  % Compute Y-axis rotation matrix (rotation by b(2)).
  c = cos(b(2)); % Cosine of Y-axis angle
  s = sin(b(2)); % Sine of Y-axis angle
  C = [c 0 -s; 0 1 0; s 0 c]; % Y-axis rotation matrix

  % Compute X-axis rotation matrix (rotation by b(1)).
  c = cos(b(1)); % Cosine of X-axis angle
  s = sin(b(1)); % Sine of X-axis angle
  B = [1 0 0; 0 c s; 0 -s c]; % X-axis rotation matrix

  % Compute combined rotation matrix (321 sequence: Z, then Y, then X).
  a = (B*C*D)'; % Transpose of B*C*D to get final rotation matrix
else
  % Throw error if input vector is not of length 3.
  error(' *** Error in ROTXYZ:  Input vector must be of length 3.');
end

return % Exit function