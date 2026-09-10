from Scheme import *
from multiprocessing import cpu_count
from time import time_ns

class Fetuccini:
    def __init__(self, reacs_in_full, reacs_from_file, inp_file='1D_PAC/hd2094_varied.inp', width=1.,
                 vary_zero=False, n_samples=10, pop_str='alpha_reacs', out_path='1D_PAC/fetuccini', command='./',
                 max_processes=12, silent=False, start_id=0, params_to_vary=['alpha'], run_rn=True,
                 mode='lognormal', use_all_reacs=True, max_width=3., neg_variation=False,
                 include_changes=True):
        self.max_width, self._width = max_width, width
        self._include_changes = include_changes
        self.out_path = self._handle_dir(out_path)
        self.reacs_files = self.cook_fetuccini(reacs_in_full, reacs_from_file, inp_template=inp_file,
                                         vary_zero=vary_zero, n_samples=n_samples,
                                          pop_str=pop_str, start_id=start_id, neg_variation=neg_variation,
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
                 use_all_reacs=False, not_include=[],
                 stop_variation=[]):
        aux_spaces = '                                                    '
        max_width = self.max_width
        width = self._width
        if mode == 'lognormal':
            mu = 0
            # Seed the random number generator with the current time in nanoseconds
            seed = int(time_ns())
            rng = np.random.default_rng(seed)
            samples = rng.normal(mu, self._width, int(10 * len(reacs_in_full)))
            filtered_gauss = samples[(samples >= mu - self.max_width) & (
                                        samples <= mu + self.max_width)]
            filtered_gauss = rng.permutation(filtered_gauss)[:len(reacs_in_full)]

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
                        new_reacs_df.at[ind, param+'_change'] = log_normal

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

            this_line = reacs_from_file[id_]
            aux_line= this_line.split('\n')
            if not i in stop_variation:
                line_i = aux_line[0].split(':')
                start_i = line_i[0][:-low_lim_int]
                end_i = ' '
                rtype = int(new_reacs_df.at[i, 'rtype'])
                new_alpha, new_beta, new_gamma = "{:.2E}".format(new_reacs_df.at[i, 'alpha']
                                                    ), "{:.2f}".format(new_reacs_df.at[i, 'beta']
                                                    ), "{:.0f}".format(new_reacs_df.at[i, 'gamma'])
            
                string_rtype = format_line(rtype, total_length=4)
                to_change = format_line(new_alpha+'     '+new_beta, new_gamma, total_length=26)

                if abs(rtype) in [800, 811, 812, 900, 911]:
                    new_line = ''
                    second_line = lines[1]
                    alpha_inf, beta_inf, gamma_inf = second_line.split()[:3]
                    new_alpha_inf = "{:.2E}".format(new_reacs_df.at[i, 'alpha_change']*float(alpha_inf))
                    new_line += aux_spaces
                    new_line += format_line(new_alpha_inf+'     '+beta_inf, gamma_inf, 
                                total_length=26) + aux_spaces[:31]+': k_infinite'
                    lines[1] = new_line

                elif 30 < abs(rtype) and abs(rtype) < 45:
                    for i, line in zip(range(len(lines)-2), lines[1:-1]):
                        alpha_plog, beta_plog, gamma_plog = line.split()[:3]
                        pref = line.split()[-1]
                        new_alpha_plog = "{:.2E}".format(new_reacs_df.at[i, 'alpha_change']*float(alpha_plog))
                        new_line = aux_spaces
                        new_line += format_line(new_alpha_plog+'     '+beta_plog, gamma_plog, 
                                total_length=26) + aux_spaces[:31]
                        new_line +=': '+pref
                        lines[i+1] = new_line

                #else: abs(rtype) in [3, 20, 802, 902]

                for st in line_i[1:]:
                    end_i+=':'+st
                for line in lines[1:]:
                    end_i+='\n'+line

                if not self._include_changes:
                    new_line = start_i +'    '+ to_change+' '+ string_rtype+ end_i
                else:
                    str_change = "{:.2E}".format(new_reacs_df.at[i, 'alpha_change'])
                    new_line = start_i +'    '+ to_change+' '+ str_change+ end_i
            else:
                new_line = this_line
            new_reacs_dict[id_] = new_line
            varied_file.write(new_line)
        varied_file.close()
        del new_reacs_df
        return new_reacs_dict
    
    def cook_fetuccini(self, reacs_in_full, reacs_from_file, pop_str='varied_reacs', inp_template = '1D_PAC/hd2094_varied.inp', 
                        vary_zero=False, params_to_vary = ['alpha'], n_samples=5, start_id=0, 
                        mode='lognormal', use_all_reacs=False, neg_variation=False):
        out_path = self.out_path
        templates = [open(inp_template, 'r').readlines()]
        reacs_files = []
        if neg_variation:
            fix_reacs = [id_ for id_ in range(len(reacs_in_full))]
            print (fix_reacs)
            fix_reacs = np.repeat(fix_reacs, n_samples)
            print (fix_reacs)
            
            n_samples *= len(fix_reacs)
        for j in range(n_samples):
            fname = out_path / Path(f'{pop_str}_{j}.reac')
            aux_sch = Scheme([], [], [], search_first=False, planets=['hd2094'], scheme_id=j+start_id,
                            out_path=out_path, pop_str=pop_str)
            aux_sch._write_scheme(fname, change_spec_file=False, templates=templates)
            new_reacs_file = self.vary_params(reacs_in_full, reacs_from_file, fname=fname,
                                              vary_zero=vary_zero, params_to_vary=params_to_vary,
                                              mode=mode, use_all_reacs=use_all_reacs, stop_variation=[] if not neg_variation else fix_reacs[i])
            reacs_files.append(new_reacs_file)
        return reacs_files
    


if __name__ == "__main__":
    num_cores = cpu_count()
    n_samples = int(10*num_cores)
    n_samples = 2
    strtime = datetime.now().strftime("var-%Y:%m:%d:%H:%M-")
    initialize(network='venot', fetuccini_mode=True)
    from darwen_functions import reacs_in_full, reacs_from_file
    Fetuccini(reacs_in_full, reacs_from_file, n_samples=n_samples, run_rn=True, max_processes=num_cores,
              pop_str=strtime, out_path = '1D_PAC/fetuccini_oct', inp_file='1D_PAC/hd2094_varied.inp',)
    
