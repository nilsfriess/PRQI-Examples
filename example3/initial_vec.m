function vec = initial_vec(x, n_osc, R)
  assert(n_osc > 0);

  value  = 1;
  period = (2 * R) / n_osc;

  vec = -value + 2*value*(mod(x - period/2, period) < 0.5 * period);
  vec = vec .* (x < R) .* (x > 0.1);
  vec = vec(:);
end
