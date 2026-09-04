function v = fm_env (key, v, parse)
% FM_ENV  Override v with the environment variable FM_<key> when it is set, converted by the optional handle
  s = getenv (['FM_' key]);
  if (isempty (s)); return; end
  if (nargin > 2); v = parse (s); else; v = s; end
end
