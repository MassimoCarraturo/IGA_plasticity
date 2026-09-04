function write_tables (R, method_tags, max_levels, results_dir, bench)
%WRITE_TABLES  Console + .dat output for the front-kink metric studies.
%
%   write_tables (R, method_tags, max_levels, results_dir, bench)
%
% R is the struct of results keyed '<tag>_lev<N>' with the fields produced by
% front_kink_metrics plus ndof/nel/wall.  Also reports whether all operators
% shared the same mesh, which is the premise of the comparison.

  for il = 1:numel(max_levels)
      lev = max_levels(il);
      fprintf('\n\n===== %s : FRONT METRICS, max_level = %d =====\n', upper(bench), lev);
      fprintf('%-10s %6s %10s %10s %8s %10s %8s %10s %8s %8s\n', ...
          'Operator','DOFs','E_glob','E_band','E_inf%','E_slope','J','V_yield','TV','wall');
      fid = fopen(fullfile(results_dir, sprintf('front_metrics_%s_lev%d.dat', bench, lev)), 'w');
      fprintf(fid, 'method\tndof\tnel\tE_glob\tE_band\tE_inf\tE_slope\tJ\tV_yield\tTV\twall\n');
      nels = nan(1, numel(method_tags));
      for im = 1:numel(method_tags)
          key = sprintf('%s_lev%d', method_tags{im}, lev);
          if ~isfield(R, key); continue; end
          t = R.(key);  nels(im) = t.nel;
          fprintf('%-10s %6d %10.3e %10.3e %8.2f %10.3e %8.3f %10.3e %8.3f %8.1f\n', ...
              method_tags{im}, t.ndof, t.E_glob, t.E_band, t.E_inf, t.E_slope, ...
              t.J, t.V_yield, t.TV, t.wall);
          fprintf(fid, '%s\t%d\t%d\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.8e\t%.3f\n', ...
              method_tags{im}, t.ndof, t.nel, t.E_glob, t.E_band, t.E_inf, ...
              t.E_slope, t.J, t.V_yield, t.TV, t.wall);
      end
      fclose(fid);
      u = unique(nels(~isnan(nels)));
      if numel(u) == 1
          fprintf('  [mesh check] identical mesh for all operators (nel = %d)\n', u);
      else
          fprintf('  [mesh check] WARNING meshes differ: nel = %s\n', mat2str(nels));
      end
  end
end
