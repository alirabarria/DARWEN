from Scheme import *

class Fetuccini:
    def __init__(self, reacs_in_full, reacs_from_file, inp_file='1D_PAC/hd2094_varied.inp', width=1.,
                 vary_zero=False, n_samples=10, pop_str='alpha_reacs', out_path='1D_PAC/fetuccini', command='./',
                 max_processes=12, silent=False, start_id=0, params_to_vary=['alpha'], run_rn=True,
                 mode='lognormal', use_all_reacs=True, max_width=3.):
        self.max_width, self._width = max_width, width
        self.out_path = self._handle_dir(out_path)
        self.reacs_files = self.cook_fetuccini(reacs_in_full, reacs_from_file, inp_template=inp_file,
                                         vary_zero=vary_zero, n_samples=n_samples,
                                          pop_str=pop_str, start_id=start_id,
                                          params_to_vary=params_to_vary, mode=mode, use_all_reacs=use_all_reacs)
        if run_rn:
            run_pac(pop_str=pop_str, silent=silent, max_processes=max_processes, command=command,
                hide_output=False, out_path=out_path, out_file=None)
        

    def _handle_dir(self, direc, abs = False):
        direc = Path(direc) if not isinstance(direc, Path) else direc
        if not abs:
            direc.resolve()
        else:
            direc.absolute()
        direc.mkdir(parents=True, exist_ok=True)
        return direc
        


    def vary_params(self, reacs_in_full, reacs_from_file, vary_zero=False, 
                    params_to_vary = ['alpha'], fname='1D_PAC/varied_reacs.reac', 
                 low_lim_int=36, mode='lognormal',
                 use_all_reacs=False, not_include=[]):
        max_width = self.max_width
        width = self._width
        if mode == 'lognormal':
            mu = 0
            samples = np.random.normal(mu, self._width, int(10*len(reacs_in_full)))
            filtered_gauss = samples[(samples >= mu - self.max_width) & (
                                        samples <= mu + self.max_width)]
            filtered_gauss = np.random.permutation(filtered_gauss)[:len(reacs_in_full)]

        new_reacs_df = reacs_in_full.copy()
        new_reacs_dict = {} 
        for ind, row in reacs_in_full.iterrows():
            for param in params_to_vary:
                if row[param] != 0 or vary_zero:   
                    if mode == 'lognormal':
                        N = filtered_gauss[ind]
                        F = row['unc_max']
                        log_normal = np.exp(np.log(F)*N)
                        new_param = log_normal*row[param]

                    else:
                        sign = np.sign(row[param])
                        sample_gauss = np.random.normal(0, width)
                        new_param = (1+sample_gauss)*row[param]
                        if sign != np.sign(new_param) and sign != 0:
                            new_param = sign*np.abs(new_param)
                    new_reacs_df.at[ind, param] = new_param

        varied_file = open(fname, 'w')
        for i in range(len(reacs_in_full)):
            if not use_all_reacs:
                id_ = reacs_in_full.at[i, 'id']
            else:
                id_ = i

            if id_ in not_include:
                continue

            line_i = reacs_from_file[id_].split(':')
            start_i = line_i[0][:-low_lim_int]
            end_i = ' '
            for st in line_i[1:]:
                end_i+=':'+st
            rtype = int(new_reacs_df.at[i, 'rtype'])
            new_alpha, new_beta, new_gamma = "{:.2E}".format(new_reacs_df.at[i, 'alpha']
                                                ), "{:.2f}".format(new_reacs_df.at[i, 'beta']
                                                ), "{:.0f}".format(new_reacs_df.at[i, 'gamma'])
            string_rtype = format_line(rtype, total_length=4)
            to_change = format_line(new_alpha+'     '+new_beta, new_gamma, total_length=26)
            new_line = start_i +'    '+ to_change+' '+ string_rtype+ end_i
            new_reacs_dict[id_] = new_line
            varied_file.write(new_line)
        varied_file.close()
        del new_reacs_df
        return new_reacs_dict
    
    def cook_fetuccini(self, reacs_in_full, reacs_from_file, pop_str='varied_reacs', inp_template = '1D_PAC/hd2094_varied.inp', 
                        vary_zero=False, params_to_vary = ['alpha'], n_samples=5, start_id=0, 
                        mode='lognormal', use_all_reacs=False):
        out_path = self.out_path
        templates = [open(inp_template, 'r').readlines()]
        reacs_files = []
        for j in range(n_samples):
            fname = out_path / Path(f'{pop_str}_{j}.reac')
            aux_sch = Scheme([], [], [], search_first=False, planets=['hd2094'], scheme_id=j+start_id,
                            out_path=out_path, pop_str=pop_str)
            aux_sch._write_scheme(fname, change_spec_file=False, templates=templates)
            new_reacs_file = self.vary_params(reacs_in_full, reacs_from_file, fname=fname,
                                              vary_zero=vary_zero, params_to_vary=params_to_vary,
                                              mode=mode, use_all_reacs=use_all_reacs)
            reacs_files.append(new_reacs_file)
        return reacs_files
    


if __name__ == "__main__":

    initialize(network='venot', fetuccini_mode=True)
    Fetuccini(reacs_in_full, reacs_from_file, n_samples=1000, run_rn=True, max_processes=36,
              pop_str='venot', out_path = '1D_PAC/fetuccini', inp_file='1D_PAC/hd2094_varied.inp',)
    
