import matplotlib.pyplot as plt

from Scheme import *

class Darwen:
    def __init__(self, mutation_rate = 0.1, mutation_prob=.1, crossover_prob=.1, max_iter=50, silent=False, 
                temp_path='temp/', imp_molecs=imp_molecs, all_molecs=all_molecs, 
                templates=['1D_PAC/reduc_p10_03.inp'], full_dirs=['1D_PAC/conv_test'], 
                fitness_metric='max_percent', big_loss_def=1e5,
                weight_by='molecs', weight=[1e-3], push_imp_molecs=True, len_reacs=958,
                molecs_low_limit=38, molecs_high_limit=55, initial_scheme_dir='1D_PAC/opt_sch_1217',
                reacs_df=reacs_in_full, logarithmic_growth=False, n_fittest=15,
                init_population=300, initial_kill_rate=0, initial_mutation_rate=0.2, n_crossover=1e3,
                mod_elitism=5, mod_curious = 10, pop_prefix='pop', mutate_mols=True, n_population=15,
                kill_rate=0.05, crossover_tries = 2, search_init='', start_changing_mols=True,
                initial_reacs=[], initial_loss=[1e3,1e3,1e3], min_is_better=True, measure_only_imp=True,
                stop_loss=0.03, kill_rate_factor=0.9, planets=['hd2094', 'hd1897'], weight_fn=np.max,
                get_initial_loss=True, abs_threshold=1e5, 
                rosetta_path='/Users/arturolira/Desktop/phd/opt',
                boost_new_mols_factor=0.5, boost_mols=True, search_schemes_first=True,
                seed=-1, saved_schemes_path='saved_schemes/', master_df_name='master_index.csv',
                master_cols=['mols_len', 'ids_len', 'ids', 'molecs', 'loss', 'delta_loss', 'losses', 'abs_losses', 
                             'all_rel_losses', 'major_loss',
                            'directory', 'fname_key', 'parent', 'len_pool', 'date', 'loss_type', 'worst_mols', 
                             'weights', 'planets', 'photo_flag', 'single_planet_flag', 'curiosity', 'genes',
                             'model_dict'],
                save_schemes_loss_threshold=.2, ids_df_name='ids_df.csv', update_ids_df=False,
                smart_choice_mode='sometimes', max_smart_ids_per_reac=3, filter_ids_threshold=1e3,
                curious_pool_ratio=1/20, filter_ids_mode='sum', plt_show=False, hide_output=False,
                max_comb_crossover=200, photo_flag=False, single_planet_flag=False, curiosity_factor=1,
                max_curious_ids=3, always_deterministic=True, n_before_wheel=2,
                extended_range_flag=False, add_per_delta=5, change_imp_molecs_flag=False):
        
        """
        Initializes the Darwen class, which is designed to optimize chemical reaction schemes using a genetic algorithm.

        Parameters:
            mutation_rate (float): The base rate at which mutations occur during the genetic algorithm process.
            max_iter (int): Maximum number of iterations the genetic algorithm will run.
            silent (bool): If True, suppresses most of the output messages to keep the execution quiet.
            temp_path (str): Base directory path where output files will be stored.
            imp_molecs (list): List of important molecules to prioritize in the optimization process.
            all_molecs (list): List of all molecules considered in the chemical reaction schemes.
            templates (list of str): Paths to template files used as starting points for generating new chemical reaction schemes.
            full_dirs (list of str): Directories containing full optimization datasets for evaluation.
            fitness_metric (str): Metric used to evaluate the fitness of a chemical reaction scheme.
            big_loss_def (float): Default large loss value used for penalization in the optimization process.
            weight_by (str): Criterion for weighting loss calculations, typically by molecule importance.
            weight (list of float): Weights applied to the fitness calculation.
            push_imp_molecs (bool): If True, emphasizes the importance of molecules specified in `imp_molecs`.
            len_reacs (int): The total number of reactions considered in the optimization space.
            molecs_low_limit (int), molecs_high_limit (int): Range of molecule counts to consider for valid chemical schemes.
            initial_scheme_dir (str): Directory containing initial chemical scheme configurations.
            reacs_df (pd.DataFrame): DataFrame containing reaction data used in the optimization process.
            logarithmic_growth (bool): Determines if the population growth should be logarithmic.
            n_fittest (int): Number of best schemes to retain between iterations.
            init_population (int): Initial population size for the genetic algorithm.
            initial_kill_rate (float), initial_mutation_rate (float): Initial rates for killing and mutating schemes in the population.
            n_crossover (int): Number of crossover operations to perform per iteration.
            n_elitism (int): Number of elite schemes to automatically carry over to the next iteration.
            pop_prefix (str): Prefix used for naming populations during the optimization process.
            mutate_mols (bool): If True, allows mutation operations to alter molecule configurations.
            n_population (int): Number of schemes in the population for each iteration.
            kill_rate (float), crossover_tries (int): Parameters controlling the intensity of killing and crossover operations.
            search_init (str): Initial search string for finding existing schemes.
            start_changing_mols (bool): If True, starts the optimization process with molecule configuration changes enabled.
            initial_reacs (list), initial_loss (list): Initial reactions and their associated loss values for starting the optimization.
            min_is_better (bool): Indicates if a lower loss value is considered better.
            measure_only_imp (bool): If True, measures fitness based only on important molecules.
            stop_loss (float): Loss threshold for early optimization termination.
            kill_rate_factor (float): Factor adjusting the intensity of scheme culling.
            planets (list of str): Planetary conditions for scheme evaluation.
            weight_fn (function): Function to aggregate losses (np.mean, np.max, etc.).
            get_initial_loss (bool): Calculates initial loss values for starting schemes.
            abs_threshold (float): Absolute loss threshold for scheme viability.
            code_path (str): Path to the main codebase or execution environment.
            rosetta_path (str): Path to the Rosetta directory or essential data set.
            boost_new_mols_factor (float): Factor for introducing new molecules.
            boost_mols (bool): Enables active introduction of new molecules.
            search_schemes_first (bool): Prefers reusing existing schemes if possible.
            seed (int): Seed for random number generation, affecting reproducibility.
            saved_schemes_path (str): Path for saving optimized schemes.
            master_df_name (str): Name for the master DataFrame CSV file logging scheme details.
            master_cols (list of str): Column names for the master DataFrame.
            save_schemes_loss_threshold (float): Loss threshold below which schemes are saved.

        """
        self.params= {'silent':silent, 'crossover_prob': crossover_prob,
                     'templates': templates, 'full_dirs': full_dirs, 'mutation_rate': mutation_rate,
                     'fitness_metric':fitness_metric, 'big_loss_def':big_loss_def, 'weight_by':weight_by,
                     'weight':weight, 'push_imp_molecs':push_imp_molecs, 'len_reacs':len_reacs,
                     'molecs_low_limit':molecs_low_limit, 'molecs_high_limit':molecs_high_limit,
                     'logarithmic_growth':logarithmic_growth, 'n_fittest': n_fittest,
                     'initial_kill_rate': initial_kill_rate, 'init_population': init_population,
                      'initial_mutation_rate':initial_mutation_rate, 'n_crossover':n_crossover,
                      'mod_elitism':mod_elitism, 'pop_prefix': pop_prefix, 'mutate_mols':mutate_mols,
                      'n_population':n_population, 'mutation_prob': mutation_prob,
                      'kill_rate':kill_rate, 'crossover_tries':crossover_tries,
                      'search_init': search_init, 'start_changing_mols':start_changing_mols,
                      'max_iter': max_iter, 'planets': planets,
                      'stop_loss': stop_loss, 'kill_rate_factor': kill_rate_factor,
                      'weight_fn': weight_fn, 'abs_threshold': abs_threshold,
                       'boost_new_mols_factor':boost_new_mols_factor,
                      'boost_mols':boost_mols, 'search_schemes_first': search_schemes_first,
                      'master_cols':master_cols,
                      'master_df_name': master_df_name, 
                      'save_schemes_threshold': save_schemes_loss_threshold,
                      'initial_scheme_dir':initial_scheme_dir, 'smart_choice_mode':smart_choice_mode,
                      'max_smart_ids_per_reac':max_smart_ids_per_reac, 
                      'filter_ids_threshold':filter_ids_threshold,
                      'curious_pool_ratio':curious_pool_ratio, 'filter_ids_mode':filter_ids_mode,
                      'plt_show': plt_show, 'max_comb_crossover': max_comb_crossover,
                      'photo_flag':photo_flag, 'single_planet_flag':single_planet_flag,
                      'curiosity_factor': curiosity_factor, 'mod_curious':mod_curious,
                      'max_curious_ids':max_curious_ids, 'n_before_wheel': n_before_wheel,
                      'always_deterministic': always_deterministic,
                      'extended_range_flag': extended_range_flag,
                      'add_per_delta': add_per_delta, 
                      'change_imp_molecs_flag': change_imp_molecs_flag
                      
                     }
        self._hide_output = hide_output
        if hide_output:
            output_file = (Path('1D_PAC') / Path(f'output_{date.today()}.txt')).resolve()
            # Ensure the parent directory exists
            output_file.parent.mkdir(parents=True, exist_ok=True)

            # Create the file
            output_file.touch(exist_ok=True)
            self._out_file = open(output_file, 'a')
        else:
            self._out_file = None
        self._seed=seed
        self.imp_molecs = imp_molecs
        self.all_molecs = all_molecs
        self.measure_only_imp = measure_only_imp
        self.min_is_better = min_is_better
        not_imp_molecs_mask = np.array([not m in imp_molecs for m in all_molecs])
        self.not_imp_molecs = self.all_molecs[not_imp_molecs_mask]

        self.reacs_df = reacs_df
        self.best_loss, self.old_losses = ( big_loss_def*np.ones((n_fittest, 3)) for i in range(2))

        
        self.best_schemes, self.best_schemes_log, self.best_loss_log, self.optimate_schemes = ([] for i in range(4))
        self._ideal_dfs, self._templates, self.optimate_losses = ([] for i in range(3))
        
        self.initial_loss = big_loss_def
        self._run_counter = 0

        self.schemes_path, self.temp_path, self.rosetta_path = self._handle_dir(saved_schemes_path
                                                            ), self._handle_dir(temp_path
                                                            ), self._handle_dir(rosetta_path+'/'+str(date.today()))
        split_ids_name = ids_df_name.split('.')
        ids_df_name_file = split_ids_name[0]
        for pl in planets:
            ids_df_name_file+= '_'+pl

        ids_df_name = ids_df_name_file+'.'+split_ids_name[1]

        self._master_df_path, self._ids_df_path = self.schemes_path / Path(master_df_name), self.schemes_path / Path(ids_df_name)

        if self._master_df_path.exists():
            self.master_df=pd.read_csv(self._master_df_path, usecols=master_cols)
        else:
            self.master_df=pd.DataFrame(columns=master_cols)

        if self._ids_df_path.exists() and not update_ids_df:
            self.ids_df=pd.read_csv(self._ids_df_path)
        elif update_ids_df:
            self.construct_ids_df()
        else:
            self.ids_df = []
        
        
        assert len(self.params['templates']) == len(self.params['full_dirs'])

        for dir_ in self.params['full_dirs']:
            ideal_df, cons_err, _err = out_to_dfs(dir_=dir_, error_purposes=True)
            if cons_err:
                print('Are you sure the loss directory is ideal?', file=self._out_file)
                if self._hide_output: self._out_file.flush()
            self._ideal_dfs.append(ideal_df)
            
        for inp_template in self.params['templates']:
            self._templates.append(open(inp_template, 'r').readlines())

        
        if search_init == '' and len(initial_reacs)==0:
            self.initial_scheme = self.get_initial_scheme(initial_scheme_dir, 
                                                        path='', get_loss=get_initial_loss)
            if get_initial_loss:
                self.initial_loss = self.initial_scheme.loss
            
        elif len(initial_reacs)>0:
            self.initial_scheme = Scheme(initial_reacs, self._molecs_from_ids(initial_reacs), self.master_df, search_first=self.params['search_schemes_first'],
                                        schemes_path=self.schemes_path, len_reacs=self.params['len_reacs'])
            self.initial_loss= initial_loss
            
        else:
            if self._seed > -1:
                np.random.seed(self._seed+1)
            initial_reacs = np.random.permutation(range(len_reacs))[:100]
            self.initial_scheme = Scheme(initial_reacs, self._molecs_from_ids(initial_reacs), self.master_df, search_first=self.params['search_schemes_first'],
                                        schemes_path=self.schemes_path, len_reacs=self.params['len_reacs'],)
        
        if self.min_is_better:
            self.optimate_loss = np.min([big_loss_def, self.initial_loss])
            self.optimate_losses.append(self.initial_loss)
            self.optimate_schemes.append(self.initial_scheme)
        else:
            self.optimate_loss = np.max([big_loss_def, self.initial_loss])
        
        self.optimate_losses = np.array(self.optimate_losses)
            
        
        
        if not silent:
            if get_initial_loss:
                self.print_loss(self.initial_scheme, print_rel_all=True)
            else:
                 print ('Initial loss:', self.initial_loss, file=self._out_file)
                 if self._hide_output: self._out_file.flush()

        return
    
    def construct_ids_df(self, loss_threshold=7, all_mols_threshold=0, verbose=True, planet_filter=[]):
        if not self.params['silent']:
            print ('Filling IDs DataFrame...', file=self._out_file)
            if self._hide_output: self._out_file.flush()
        self.schemes_dict = {}
        self.best_schemes_df = self.master_df[self.master_df.loss < loss_threshold]
        if len(self.best_schemes_df) ==0:
            print ('Maybe increase loss threshold?')
        if len(planet_filter) > 0:
            mask_planets = [np.any([aux_planet in aux_df for aux_planet in planet_filter]) for aux_df in self.best_schemes_df['planets']]
            self.best_schemes_df = self.best_schemes_df[mask_planets]
        
        if all_mols_threshold > 0:
            mask = [np.max(from_df_to_array(aux)) < all_mols_threshold for aux in self.best_schemes_df['all_rel_losses']]
            self.best_schemes_df = self.best_schemes_df[mask]
        if len(self.best_schemes_df) > 0:
            best_ids = self.best_schemes_df['ids']
            num_reacs = self.params['len_reacs']
            self.ids_df = pd.DataFrame(np.zeros((num_reacs, num_reacs)), columns=range(num_reacs), index=range(num_reacs))
            for i in range(num_reacs):
                if not self.params['silent'] and verbose:
                    if i%100 == 0:
                        print (int(i/num_reacs*100), '%', end=' ', file=self._out_file)
                        if self._hide_output: self._out_file.flush()
                aux_num = str(float(i))
                mask_id = [(' '+aux_num in aux or '['+aux_num in aux) for aux in best_ids]
                this_i_ids = best_ids[mask_id]
                self.schemes_dict[i]=this_i_ids
                for j in range(num_reacs):
                    if i < j:
                        continue
                        #self.ids_df.at[i, j] = self.ids_df.at[j,i]
                    elif i==j:
                        self.ids_df.at[i,i] = len(this_i_ids)
                    else:
                        aux_num = str(float(j))
                        self.ids_df.at[i,j] = np.sum([(' '+aux_num in aux or '['+aux_num in aux) for aux in this_i_ids])
            self.ids_df.to_csv(self._ids_df_path, index=False)
            print ('100 %', file=self._out_file)
        else:
            self.ids_df = []
            print ('Failed to build ids df for planets:', planet_filter, file=self._out_file)
        if self._hide_output: self._out_file.flush()

                
    
    def _handle_dir(self, direc, abs = False):
        direc = Path(direc) if not isinstance(direc, Path) else direc
        if not abs:
            direc.resolve()
        else:
            direc.absolute()
        direc.mkdir(parents=True, exist_ok=True)
        return direc


    def clean_dir(self, directory='', str_to_remove=['pop', 'output_pop']):
        
        """
        Removes files or directories within the given directory that start with any of the specified strings.

        Parameters:
        - directory (Path): The directory to clean, specified as a Path object.
        - str_to_remove (list): A list of strings. Files or directories starting with these strings will be removed.
        """
        if directory=='':
            directory = self.temp_path
        # Iterate over files and directories in the specified directory
        for item in directory.iterdir():
            # Check if the item's name matches the pattern to be removed
            if any(item.name.startswith(prefix) for prefix in str_to_remove):
                # Remove the file or directory
                if item.is_file():
                    item.unlink()  # Removes a file
                elif item.is_dir():
                    shutil.rmtree(item)  # Removes a directory and all its contents
                    
    def print_loss(self, scheme, round_=6, print_rel_all = True, shout=False):

        if self.params['silent'] and not shout:
            return
        try:
            print (file=self._out_file)
            print('.........................................................................................................', file=self._out_file)
            if not print_rel_all:
                print ('Planet','       ', self.params['fitness_metric'],'        ', 'Absolute error','       ', 'Worst molec', file=self._out_file)
            else:
                print ('Planet','       ', self.params['fitness_metric'],'        ', 'Absolute error','       ', 'Worst molec',
                      '      ', self.params['fitness_metric']+' for all molecs', '    ...for major species', file=self._out_file)
            print('.........................................................................................................', file=self._out_file)
            if self._hide_output: self._out_file.flush()
            for id_planet in range(len(self.params['planets'])):
                if not print_rel_all:
                    print (self.params['planets'][id_planet],'  ', scheme.losses[id_planet],'  ', 
                       scheme.max_err[id_planet],'       ', scheme._worst_mols[id_planet], file=self._out_file)
                else:
                    print (self.params['planets'][id_planet],'  ', scheme.losses[id_planet],'  ', 
                       scheme.max_err[id_planet],'       ', scheme._worst_mols[id_planet], '       ', 
                           scheme._rel_all[id_planet], scheme._major_loss[id_planet], file=self._out_file)

            print('.........................................................................................................', file=self._out_file)
            move='                                                                                   '
            print (move, 'Number of molecs:', len(scheme.molecs), file=self._out_file)
            print (move[:-4],'Imp molecs present?', np.all([m in scheme.molecs for m in self.imp_molecs]), file=self._out_file)
            if self.params['weight_by']=='molecs':
                print (move[:-65], 'Total loss:', scheme.loss ,'=', self.params['weight_fn'](scheme.losses), '+', 
                       self.params['weight'][0],'*', len(scheme.molecs), file=self._out_file)
            elif self.params['weight_by']=='molecs&allabs':
                print (move[len(move)-65+3*round_+(round_-1):],'Total loss:', round (self.params['weight_fn'](scheme.losses), round_), '+', 
                       round(self.params['weight'][0] * len(scheme.molecs), round_), '+',
                      round(self.params['weight'][1] * self.params['weight_fn'](scheme.max_err), round_),'=', scheme.loss, file=self._out_file)
                
            elif self.params['weight_by'] == 'molecs&relmajor':
                print (move[len(move)-65+3*round_+(round_-1):],'Total loss:', round (self.params['weight_fn'](scheme.losses), round_), '+', 
                       round(self.params['weight'][0] * len(scheme.molecs), round_), '+',
                      round(self.params['weight'][1] * self.params['weight_fn'](scheme._major_loss), round_),'=', scheme.loss, file=self._out_file)
            
            elif self.params['weight_by'] == 'molecs&logmajor':
                print (move[len(move)-65+3*round_+(round_-1):], 'Total loss:', round (self.params['weight_fn'](scheme.losses), round_), '+', 
                       round(self.params['weight'][0] * len(scheme.molecs), round_), '+',
                      round(self.params['weight'][1] * np.log10(self.params['weight_fn'](scheme._major_loss)), round_),'=', scheme.loss, file=self._out_file)

            else:
                print ('Total loss:', scheme.loss, ' weighting as', self.params['weight_by'], file=self._out_file)
            print('.........................................................................................................', file=self._out_file)
            print(file=self._out_file)
            if self._hide_output: self._out_file.flush()
        except Exception:
            print ("Couldn`t print", file=self._out_file)
            if self._hide_output: self._out_file.flush()
        return
    
    def schemes_heatmap(self, loss_thresholds = [.5, 1, 2], apply_context=True, n_pops=13, init_str='',
                        planets=[], fname='genes.pdf', filter_by_save=False, use_real_losses=False,
                        colors=['royalblue', 'navy', 'blue', 'lightblue'], cmap_color='Oranges',
                        legend_coords=(1., 1.1), adjust = [0.1, 0.92, 0.2, 0.85],
                        gridspec_kw={"height_ratios":[.1, .1, 1], "width_ratios": [1,.05],
                                    "wspace":0.05, "hspace":0.05},
                                    filter_by_date=True):
        plt.close('all')
        import seaborn as sns
        if apply_context:   
            sns.set_context('talk')
        if filter_by_save:
            useful_df = self.master_df[self.master_df != 'Not saved']
        else:
            useful_df = self.master_df

        pop_df = useful_df[[init_str in aux for aux in useful_df['fname_key']]]
        if filter_by_date:
            today = str(date.today())
            yest = get_yesterday(today)
            print (today, yest)
            pop_df = pop_df[[aux == today or aux == yest for aux in pop_df['date']]]

        for index, row in pop_df.iterrows():
            genes = row['genes']
            if isinstance(genes, str):
                genes=genes.replace(' ',',')

            try:
                better_genes = from_df_to_array(genes)
                better_genes = np.array(better_genes).astype(int)
                pop_df.at[index, 'genes'] = better_genes
            except Exception as e:
                print ('problem', index, genes)
                print(repr(e))


            
        
        
        genes_array = np.zeros((len(pop_df), 958))
        for i in range(len(pop_df['genes'].values)):
            genes_array[i] = pop_df['genes'].values[i]
        genes_df = pd.DataFrame(genes_array, columns=range(958))
        genes = []
        losses = []
        for i in range(n_pops):
            pop_str = '-'+str(i)+'_'
            mask_pop = [pop_str in aux for aux in pop_df['fname_key']]
            this_pop_genes = genes_df[mask_pop]
            this_pop_df = pop_df[mask_pop]
            if len(this_pop_df) == 0:
                break
            genes_av = this_pop_genes.mean().values

            if use_real_losses:
                aux_losses = []
                for larray in this_pop_df['losses']:
                    aux_losses.append(from_df_to_array(larray)[0])
                aux_losses = np.array(aux_losses)
                count_min = np.sum(aux_losses<loss_thresholds[0] )/len(aux_losses)
                if len(aux_losses) > 0:
                    min_loss = np.min(aux_losses)
                else:
                    min_loss = self.params['big_loss_def']
                
                losses.append([min_loss, count_min])
            else:
                aux_add = []
                aux_add.append(this_pop_df['loss'].min())
                for lt in loss_thresholds:
                    aux_add.append(np.sum(this_pop_df['loss']<lt)/len(this_pop_df))
                losses.append(aux_add)
            genes.append(genes_av)
            
            

        genes = np.array(genes)
        losses = np.array(losses)
        fig, axes = plt.subplots(nrows=3, ncols=2, figsize=(7,13),  gridspec_kw=gridspec_kw)
        ax_loss = axes[0,0]
        ax_count = axes[1,0]
        ax = axes[2,0]
        cax = axes[2,1]
        
        axes[0,1].set_visible(False)
        axes[1,1].set_visible(False)
        trans_genes = genes.transpose()
        
        ax_loss.plot(range(len(losses)), losses[:, 0], '*-', color=colors[0])

        for l in range(len(loss_thresholds)):
            lt = loss_thresholds[l]
            ax_count.plot(range(len(losses)), losses[:, l+1]*100, '.-', color=colors[l+1], label=lt)
        ax_count.legend(loc='upper left', bbox_to_anchor=legend_coords)
        ax_count.set_ylabel(r'$\% <$')
        ax_loss.set_ylabel(r'$\varphi$')

        sns.heatmap(trans_genes, cmap=cmap_color, ax=ax, cbar=False)
        fig.colorbar(ax.get_children()[0], cax=cax, orientation="vertical")
        ax.set_xlabel('Generation ID (~convergence time)' )
        ax.set_ylabel('Reaction ID')
        if init_str=='':
            to_title = 'Genes evolution for '+str(planets)
        else:
            to_title = 'Genes evolution for a single run of '+str(planets)
        ax_loss.set_title(to_title)
        plt.setp(ax_loss.get_xticklabels(), visible=False)
        if len(adjust) > 0:
            plt.subplots_adjust(bottom=adjust[0], top=adjust[1], left=adjust[2], right=adjust[3])
        fig.savefig(fname)
        if self.params['plt_show']:
            plt.show()
            
    def best_schemes_hist(self, scheme_list, figsize=(13,7), merge_fn=np.max, min_len=30, max_len=55, big_loss_def=1e5, 
                          fname='1D_PAC/best_hist_', 
                    plot_threshold=-1, colors=['coral', 'gray', 'darkgreen', 'green']):
        if plot_threshold == -1:
            plot_threshold = big_loss_def
        plt.close('all')
        fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=figsize, sharex=True)
        best_per_mol = {}
        for l in range(min_len, max_len+1):
            best_per_mol[l] = [big_loss_def, 0, 0, 0]
        emp_min_len = 1e3
        emp_max_len = 0
        for scheme in scheme_list:
            mols = scheme.molecs
            ids = scheme.ids
            new_len_pool = scheme._len_ids_pool + len(ids)
            real_loss = merge_fn(scheme.losses)
            if len(mols) < emp_min_len:
                emp_min_len = len(mols)
            elif len(mols) > emp_max_len:
                emp_max_len = len(mols)
            if len(mols) in best_per_mol.keys():
                best_loss, count_ids, count_schemes, len_pool = best_per_mol[len(mols)]
            else:
                best_per_mol[len(mols)] = [big_loss_def, 0, 0, 0]
            if real_loss < best_loss:
                best_loss = real_loss
                count_ids = len(ids)
                len_pool = new_len_pool
            count_schemes+=1
            best_per_mol[len(mols)] = [best_loss, count_ids, count_schemes, len_pool]
        for lm in range(int(emp_min_len), int(emp_max_len+1)):
            if lm in best_per_mol.keys():
                best_loss, count_ids, count_schemes, len_pool = best_per_mol[lm]
            else:
                best_per_mol[lm] = [big_loss_def, 0, 0, 0]

            if count_schemes !=0 and best_loss < plot_threshold:
                ax1.plot(lm, best_loss, '*', color=colors[0])
                ax2.bar(lm, count_schemes, color=colors[1])
                ax3.bar(lm, len_pool, edgecolor=colors[2], color='None')
                ax3.bar(lm, count_ids, color=colors[3])
        ax1.set_ylabel('Best loss')
        ax2.set_ylabel('# schemes')
        ax3.set_ylabel('# best ids')
        ax3.set_xlabel('# molecs')
        if self.params['plt_show']:
            plt.show()
        fig.savefig(fname+str(date.today())+'.png')
        return best_per_mol
        

        
    def get_initial_scheme(self, dir_to_check, path='', get_loss=False, error_verbose=False):
        ids_file = list_archives(path+dir_to_check, format='.int')
        if len(ids_file)>0:
            ids_scheme = np.genfromtxt(ids_file[0])
            mols_to_return = self._molecs_from_ids(ids_scheme)
        else:
            print (path+dir_to_check, ids_file, file=self._out_file)
            if self._hide_output: self._out_file.flush()
            raise FileNotFoundError("IDS file didn't found")

        initial_scheme = Scheme(ids_scheme, mols_to_return, self.master_df, scheme_id=0, out_path=self.temp_path, 
                                search_first=self.params['search_schemes_first'], len_reacs=self.params['len_reacs'],
                                schemes_path=self.schemes_path, parent_key='big_bang', planets=self.params['planets'],
                                )
        if get_loss:
            
            init_str = dir_to_check.split('/')[1]
            initial_scheme._write_scheme(write_ids=False, check_dir=True, init_str=init_str, imp_mols_flag=self.measure_only_imp,
                                                   templates=self._templates, out_path=self.temp_path)

            for planet_dir in self.params['planets']:
                files = list_archives(self.temp_path, format=init_str, key=planet_dir+'_out')
                if len(files)==0:
                    self._run_pac(init_str)
                    break
     
            #value = <value_if_true> if <expression> else <value_if_false>
            molecs_to_measure = self.imp_molecs if self.measure_only_imp else []
            initial_scheme.get_loss(self._ideal_dfs, big_loss_def=self.params['big_loss_def'], weight_by=self.params['weight_by'], 
                                                imp_molecs=molecs_to_measure, mode=self.params['fitness_metric'], 
                                                 merge_fn=self.params['weight_fn'], #min_is_better=self.min_is_better, deprecated
                                                weight_factor=self.params['weight'], path=self.temp_path,
                                   abs_threshold= self.params['abs_threshold'], error_verbose=error_verbose)
            
        
        return initial_scheme

        
    def _molecs_from_ids(self, reacs_array):
        molecs_to_return = np.array([])
        for j in reacs_array:
            molecs_to_return = np.append(molecs_to_return, self.reacs_df.at[j, 'molecs'])
        return np.unique(molecs_to_return)
        
    def _crossover(self, n_times=1, pop_str='cpop', max_comb=100, search_schemes_first=False):
        """
        Performs crossover operations on the best schemes to generate new schemes by combining features
        of two parent schemes to create offspring with characteristics of both.

        Parameters:
        - n_times (int): Number of times to perform crossover for each pair of best schemes.
        - pop_str (str): String identifier for the population, used in naming schemes.

        Returns:
        - new_schemes (list): List of newly created schemes through crossover.
        """
        new_schemes = [] 
        c=0
        self._run_counter=0
        min_len = np.min([max_comb*n_times, n_times*len(self.best_schemes)**2])
        for _ in range(int(min_len)):
            ind0 = np.random.randint(0, len(self.best_schemes))
            ind1 = np.random.randint(0, len(self.best_schemes))
            parent0, parent1 = self.best_schemes[ind0], self.best_schemes[-ind1]
            parent_key = parent0._fname_key+'+'+parent1._fname_key
            parent_loss = np.max([parent0.loss, parent1.loss])
            for _ in range(n_times):
                # Simplified crossover logic
                crossover_point = np.random.randint(1, min(len(parent0.ids), len(parent1.ids)))
                offspring_ids = np.concatenate((parent0.ids[:crossover_point], parent1.ids[crossover_point:]))
                offspring_ids = offspring_ids  # Ensure unique IDs

                # Assuming _molecs_from_ids and validation logic are simple enough to inline here
                offspring_molecs = self._molecs_from_ids(offspring_ids)
                if self.params['molecs_low_limit'] <= len(offspring_molecs) and len(
                    offspring_molecs) <= self.params['molecs_high_limit']:
                    # Construct new scheme directly here
                    
                    new_scheme = Scheme(offspring_ids, offspring_molecs, self.master_df, parent_key=parent_key, 
                                        scheme_id=len(self.best_schemes) + len(new_schemes),
                                        out_path=self.temp_path, pop_str=pop_str, parent_loss=parent_loss,
                                        planets=self.params['planets'], len_ids_pool=len(offspring_ids), 
                                        search_first=search_schemes_first,
                                        schemes_path=self.schemes_path, len_reacs=self.params['len_reacs'],
                                       )
                    loss_there = False
                    if new_scheme._found or not search_schemes_first:
                        if search_schemes_first:
                            loss_type = self.params['fitness_metric']+'-'+self.params['weight_by']
                            loss_there = new_scheme._found_data['directory'] != 'Not saved' or new_scheme._found_data['loss_type'] == loss_type
                        else:
                            loss_there = False

                        if not loss_there:
                            new_scheme._write_scheme(write_ids=False, check_dir=True, init_str=pop_str, imp_mols_flag=self.measure_only_imp,
                                                   templates=self._templates)
                            self._run_counter+=1
                    new_schemes.append(new_scheme)
                    c+=1
            if c > max_comb:
                break
        return new_schemes, c
    
    
    def _get_reacs_pool(self, reacs_array, molec_array, force_total_pool=False):
        reacs_pool = []
        for j in range(self.params['len_reacs']):
            if not j in reacs_array or force_total_pool:
                add_to_pool = np.all([m in molec_array for m in self.reacs_df.at[j,'molecs']])
                if add_to_pool:
                    reacs_pool.append(j)
        return reacs_pool
        
    def _add_reacs(self, reacs_array, molec_array, reacs_pool, mutation_int=10, add_seed=0,
                    smart_threshold=1e5, curious_pool_ratio=.25, max_ids_per_reac=2, mode='random',
                    max_curious_ids=5, return_curiosity=False):
        c = 0
        curiosity = 0
        if len(self.ids_df) > 0:
            all_counts = self.ids_df.sum()
            curiosity_ranking = (np.max(all_counts) - all_counts) / np.max(all_counts)

            if mode=='smart':
                filtered_df = self.ids_df[self.ids_df.apply(lambda row: np.any(row > smart_threshold), axis=1)]
                for _id in np.random.permutation(reacs_array):
                    if _id in filtered_df.index:
                        row = filtered_df.loc[_id]
                        mask = [row[i] > smart_threshold and i not in reacs_array for i in range(len(row))]
                        aux_ids = np.argsort(row[mask])[::-1]
                        if len(aux_ids) > 0:
                            if len(aux_ids) > max_ids_per_reac:
                                aux_ids = aux_ids[:max_ids_per_reac]
                            add_ids = filtered_df.index[aux_ids]
                            reacs_array = np.append(reacs_array, add_ids)
                            c+=len(aux_ids)
                        if c > mutation_int:
                            break
                mutation_int = mutation_int - c

            elif mode=='curious' and max_curious_ids > 0:
                curious_pool = np.argsort(-curiosity_ranking)[:int(curious_pool_ratio*len(curiosity_ranking))]
                curious_int = np.random.randint(1, max_curious_ids)
                curious_reacs = np.random.choice(curious_pool, size=curious_int, replace=False)
                reacs_array = np.append(reacs_array, curious_reacs)
                mutation_int = mutation_int - len(curious_reacs)

        if len(reacs_pool) > mutation_int and mutation_int > 0 : #mode='random'
            if mutation_int > 0:
                adding_reacs = np.random.choice(reacs_pool, size=mutation_int, replace=False)
                reacs_array = np.sort(np.append(reacs_array, adding_reacs))
        molec_array = self._molecs_from_ids(reacs_array)
        if return_curiosity:
            if len(self.ids_df) > 0 and self.params['curiosity_factor'] != 0:
                int_array = [ int(i) for i in reacs_array]
                curiosity = self.params['curiosity_factor']*np.max(curiosity_ranking[int_array])
            return reacs_array, molec_array, curiosity
        else:
            return reacs_array, molec_array
    
    def _mutate_mols(self, reacs_array, molec_array, ensure_uniqueness=False, boost=False, boost_factor=0.5):

        if ensure_uniqueness:
            reacs_array = np.unique(reacs_array)
            molec_array= np.unique(molec_array)
        if self.params['push_imp_molecs']:
            not_imp_mask = np.array([m in self.not_imp_molecs for m in molec_array])
            molec_array = molec_array[not_imp_mask]
            correct_len = len(self.imp_molecs)
            molecs_pool = self.not_imp_molecs
        else:
            molecs_pool = self.all_molecs
            correct_len = 0
        molecs_low_limit = self.params['molecs_low_limit'] - correct_len
        molecs_high_limit = self.params['molecs_high_limit'] - correct_len
        
        if len(molec_array) > molecs_high_limit:
            new_len = molecs_high_limit - np.random.randint(1, molecs_high_limit - molecs_low_limit)
        elif len(molec_array) < molecs_low_limit:
            new_len = molecs_low_limit
        elif boost:
            mid = int((molecs_high_limit - molecs_low_limit)/2)
            if self._seed > -1:
                np.random.seed(self._seed+4)
            new_len = np.random.randint(molecs_low_limit + mid, molecs_high_limit)
        else:
            if self._seed > -1:
                np.random.seed(self._seed)
            new_len = np.random.randint(molecs_low_limit, molecs_high_limit)
        
        add_len = new_len - len(molec_array) - 1
        
        if add_len < 0 and molecs_low_limit < len(molec_array):
            if self._seed > -1:
                np.random.seed(self._seed+5)
            molec_array = np.random.choice(molec_array, size = new_len, replace=False)
            reacs_pool = self._get_reacs_pool(reacs_array, molec_array, force_total_pool=True)
            mask_allow = [aux in reacs_pool for aux in reacs_array]
            reacs_array = reacs_array[mask_allow]
            
            
        elif add_len > 0 and (len(molec_array) <= molecs_high_limit and len(molecs_pool) > add_len):
            if self._seed > -1:
                np.random.seed(self._seed+6)
            add_molecs = np.random.choice(molecs_pool, size=add_len, replace=False)
            molec_array = np.append(molec_array, add_molecs)
            if boost:
                if self.params['push_imp_molecs']:
                    aux_molecs = np.append(add_molecs, self.imp_molecs)
                else:
                    aux_molecs=add_molecs
                new_mols_pool = self._get_reacs_pool(reacs_array, aux_molecs)
                boost_int = int(boost_factor*self.params['mutation_rate']*len(new_mols_pool))
                reacs_array, molec_array = self._add_reacs(reacs_array, molec_array, new_mols_pool, mutation_int=boost_int)
        #else, molec_array stay the same
        
        if self.params['push_imp_molecs']:
            molec_array = np.append(self.imp_molecs, molec_array)
            
           
        return reacs_array, molec_array
            
            


        
        
    def mendel(self, reacs_array, molec_array, mutate_mols=True, n_population=10, mutation_rate=0, kill_rate=0, 
               verbose=False, forbid_crossover=False,
               ensure_uniqueness=False, pop_str='pop0', parent_key='', parent_loss=1e5,
               start_id=0, boost_mols=True, ensure_n_mols=True, smart_mode='sometimes', smart_threshold=1000, 
               max_threshold_ids=3, curious_pool_ratio=0.25, always_curiosity=True, forbid_self_mutation=False,
               max_curious_ids=5, hide_warning=True):
        
        """
            Applies genetic operations to evolve a population of chemical schemes towards better fitness.

        Parameters:
            reacs_array (np.array): Array of reaction IDs to be considered in the initial scheme.
            molec_array (np.array): Array of molecule IDs involved in the initial scheme.
            mutate_mols (bool): Flag to perform mutation operations on molecules. Default is True.
            n_population (int): The number of schemes to generate in the population. Default is 10.
            mutation_rate (float): The rate at which mutations occur. Default is 0.
            kill_rate (float): The rate at which schemes are culled from the population. Default is 0.
            verbose (bool): If True, prints detailed logs about the genetic operations. Default is False.
            ensure_uniqueness (bool): Ensures that only unique reactions and molecules are considered. Default is True.
            pop_str (str): Population string identifier for naming and tracking. Default is 'pop0'.
            start_id (int): Starting ID for schemes generated in this population. Helps in tracking and uniqueness.
            boost_mols (bool): If True, increases the diversity of molecules during mutation. Default is True.
            ensure_n_mols (bool): Ensures the number of molecules stays within specified limits. Default is True.
            parent_key (str): A string identifier for the parent scheme, used in logging and tracking the evolution.

        Returns:
            schemes (list): A list of Scheme objects that represent the new generation of schemes after applying genetic operations.
        """
        if ensure_uniqueness:
            if verbose:
                print ('Uniqueness?', len(reacs_array)==len(np.unique(reacs_array)), file=self._out_file)
                if self._hide_output: self._out_file.flush()
            reacs_array = np.unique(reacs_array)
            molec_array = np.unique(molec_array)
        molec_array0 = np.copy(molec_array)
        schemes = []
        reacs_array0 = np.copy(reacs_array)
        kill_rate0 = kill_rate
        mutation_rate0 = mutation_rate
        first = True
        
        
        if verbose:
            print ('Initial scheme has', len(molec_array), 'molecs and', len(reacs_array), 'reactions.', file=self._out_file)
            if self._hide_output: self._out_file.flush()

        for i in range(n_population):      
            if mutate_mols:
                reacs_array, molec_array = self._mutate_mols(reacs_array, molec_array, boost=False, 
                                                             ensure_uniqueness=ensure_uniqueness)
            
            if i < n_population*self.params['crossover_prob'] and not forbid_crossover:
                cross_scheme = self.best_schemes[np.random.randint(0, len(self.best_schemes))]
                cross_genes = cross_scheme.genes
                break_point0 = np.random.randint(0, len(cross_genes))
        
                actual_genes = cross_scheme.ids_to_genes(reacs_array, len_genes=self.params['len_reacs'])

                new_genes = np.append(cross_genes[:break_point0], actual_genes[break_point0:])
                reacs_array = np.where(new_genes)[0]
                molec_array = self._molecs_from_ids(reacs_array)
                
            reacs_pool = self._get_reacs_pool(reacs_array, molec_array)
            len_possible_reacs = len(reacs_pool)+len(reacs_array)
      
            if mutation_rate < 1 and not self.params['logarithmic_growth']:
                mutation_int = int(mutation_rate*len_possible_reacs)
            elif mutation_rate < 1:
                mutation_int = int(mutation_rate*len(reacs_pool))
            if mutation_int > len(reacs_pool): 
                if (not self.params['silent'] and first) and not hide_warning:
                    warnings.warn(f"One or more pop with len too big ({len(molec_array)}M, {len_possible_reacs} max reacs)", UserWarning)
                    first = False
                mutation_int = int(0.5*mutation_int)
                kill_rate = self.params['kill_rate_factor']*mutation_rate
                if boost_mols:
                    reacs_array, molec_array = self._mutate_mols(reacs_array, molec_array, boost=True, ensure_uniqueness=ensure_uniqueness,
                                                                boost_factor=self.params['boost_new_mols_factor'])
                
            if kill_rate != 0:
                killing_int = int(kill_rate*len(reacs_array))
                if self._seed > -1:
                    np.random.seed(self._seed+7)
                reacs_array = np.random.choice(reacs_array, size=len(reacs_array) - killing_int, replace=False)
                if verbose and i==0:
                    print ('Killing will delete ~', killing_int, 'reactions.', file=self._out_file)
                    if self._hide_output: self._out_file.flush()
                
            if verbose and i==0:
                print ('Mutation will add ~', mutation_int, 'reactions.', file=self._out_file)
                print (file=self._out_file)
                if self._hide_output: self._out_file.flush()
            
            mode='random'
            parent_key0 = parent_key
            curiosity = 0
            if (smart_mode == 'sometimes' and (i%5)==0) or smart_mode =='always':
                mode = 'smart'
                parent_key+='-smart'
            elif smart_mode == 'curious':
                mode = 'curious'
                parent_key+='-curious'

            if len(self.ids_df) == 0:
                mode='random'
                parent_key = parent_key0
          
            
            if always_curiosity or mode == 'curious':
                reacs_array, molec_array, curiosity = self._add_reacs(reacs_array, molec_array, reacs_pool, 
                                                       mutation_int=mutation_int, add_seed=i,
                                                       mode=mode, max_curious_ids= max_curious_ids,
                                                       smart_threshold=smart_threshold, 
                                                       max_ids_per_reac=max_threshold_ids,
                                                       return_curiosity=True, 
                                                       curious_pool_ratio=curious_pool_ratio)
            else:
                reacs_array, molec_array = self._add_reacs(reacs_array, molec_array, reacs_pool, 
                                                       mutation_int=mutation_int, add_seed=i,
                                                       mode=mode, max_curious_ids= max_curious_ids,
                                                       smart_threshold=smart_threshold, 
                                                       max_ids_per_reac=max_threshold_ids)

            molec_array = molec_array
            scheme_i = Scheme(reacs_array, molec_array, self.master_df, parent_key=parent_key, scheme_id=i+start_id, out_path=self.temp_path, 
                              pop_str=pop_str, parent_loss=parent_loss, len_reacs=self.params['len_reacs'],
                              planets=self.params['planets'], len_ids_pool=len(reacs_pool), search_first=self.params['search_schemes_first'],
                              schemes_path=self.schemes_path, curiosity=curiosity)
            
            n_fittest = self.params['n_fittest']
            if np.random.randint(0, int(n_population*n_fittest)
                                 ) < self.params['mutation_prob']*n_population*n_fittest and not forbid_self_mutation:
                scheme_i.self_mutate(self.params['mutation_prob'])
                scheme_i.molecs = self._molecs_from_ids(scheme_i.ids)
                _ = scheme_i.search_scheme(self.schemes_path, self.master_df)

            lmols = len(scheme_i.molecs)
            in_mols_range = (self.params['molecs_low_limit']<= lmols and lmols <= self.params['molecs_high_limit']) or not ensure_n_mols
            c = 0
            if in_mols_range:
                loss_there = False
                if scheme_i._found:
                    loss_type = self.params['fitness_metric']+'-'+self.params['weight_by']
                    loss_there = scheme_i._found_data['directory'] != 'Not saved' or scheme_i._found_data['loss_type'] == loss_type
                if not loss_there:
                    scheme_i._write_scheme(write_ids=False, check_dir=True, init_str=pop_str, imp_mols_flag=self.measure_only_imp,
                                                   templates=self._templates)
                    c+=1
            if in_mols_range:
                schemes.append(scheme_i)
            
            molec_array = molec_array0
            reacs_array = reacs_array0
            mutation_rate = mutation_rate0
            kill_rate = kill_rate0
            if mode != 'random':
                parent_key=parent_key0
            self._run_counter+=c
            
        return schemes

    
    
    def _write_schemes(self, reset_init=True, opt_id='', init_str='opt', write_all_best=True, loss_threshold=1.5):
        print ('Writing best schemes', file=self._out_file)
        if self._hide_output: self._out_file.flush()
        ind_best_scheme = -1
        if write_all_best:
            for best_scheme in self.best_schemes_log:
                self.master_df = best_scheme.save_scheme(self.master_df, templates=self._templates, imp_mols_flag=self.measure_only_imp, save_threshold=loss_threshold,
                                                        weights=self.params['weight'], photo_flag=self.params['photo_flag'],
                                                        single_planet_flag=self.params['single_planet_flag'], model_dict=self.params)
            self.master_df.to_csv(self._master_df_path, index=False)
            
            
            
        if self.min_is_better:
            ind_best_scheme = np.argmin(self.optimate_losses)
        else:
            ind_best_scheme = np.argmax(self.optimate_losses)

        self.optimate_loss = self.optimate_losses[ind_best_scheme]
        self.optimate_scheme = self.optimate_schemes[ind_best_scheme]
        
        if reset_init:
            self.initial_scheme = self.optimate_scheme
            self.initial_loss = self.optimate_losses[ind_best_scheme]

        if not self.params['silent']:
            print ('Best losses:', self.optimate_scheme.losses, file=self._out_file)
            if self._hide_output: self._out_file.flush()
            
        self.optimate_scheme._write_scheme(write_ids=True, check_dir=True, init_str=init_str+opt_id+str(self._cwrite), 
                                           imp_mols_flag=self.measure_only_imp, out_path=self.rosetta_path, 
                                           templates=self._templates)
        
        
        self._cwrite+=1
        
    def _run_pac(self, pop_str, silent=False, max_processes=6, command='./'): #command='sh '
        
        run_pac(pop_str, silent=silent, max_processes=max_processes, command=command,
                hide_output=self._hide_output, out_path=self.temp_path, out_file=self._out_file)
    
            

        

    #Main program:

    def run(self, n_iterations=-1, opt_str='opt', mutation_rate=-1, n_write=5, initialize=True, 
            multiple_runs=True, min_iteration_write=1, plot_loss=False, figsize=(7,5), deallocate=True,
            n_print_loss_detail=3, save_all_schemes=False, print_first_loss=False, max_processes=6,
            crossover_pop=2):
        """
        Executes the genetic algorithm to optimize chemical reaction schemes.

        Parameters:
            n_iterations (int): Number of iterations for the genetic algorithm. If -1, uses `max_iter` from the initialization parameters.
            opt_str (str): A string prefix for output files indicating optimization iteration.
            mutation_rate (float): Overrides the mutation rate specified during initialization.
            n_write (int): Number of iterations between writing intermediate results to disk.
            initialize (bool): If True, initializes the optimization process; if False, continues from the current state.
            multiple_runs (bool): If True, allows for multiple optimization runs.
            min_iteration_write (int): Minimum iteration number required before writing results to disk.
            plot_loss (bool): If True, generates and displays a plot of loss over iterations.
            figsize (tuple): Size of the plot if `plot_loss` is True.
            deallocate (bool): If True, clears memory of schemes that are no longer needed.
            n_print_loss_detail (int): Frequency of detailed loss printing.
            save_all_schemes (bool): If True, saves all schemes regardless of their performance.

        Returns:
            None: This method orchestrates the genetic algorithm and outputs results but does not return a value.
        """
        print (file=self._out_file)
        print('.........................................................................................................', file=self._out_file)
        print('                                       Running DARWEN                                                ', file=self._out_file)
        print('              A survival of the fittest code, but calm down, just for chemical schemes                    ', file=self._out_file)
        print('                                     Date:', str(date.today()), file=self._out_file)
        print('.........................................................................................................', file=self._out_file)
        print(file=self._out_file)
        if self._hide_output: self._out_file.flush()
        
        if mutation_rate != -1:
            self.params['mutation_rate'] = mutation_rate
            

        if n_iterations <= 0:
            n_iterations = self.params['max_iter']
            
        if initialize:
            self.old_schemes = []
            self.best_schemes = []
            pop_str = self.params['pop_prefix']+'0'
            self._cwrite = 0
            if print_first_loss:
                self.print_loss(self.initial_scheme)
            if self.params['curious_pool_ratio'] == 0 or self.params['max_curious_ids'] == 0:
                smart_mode = 'random'
            else:
                smart_mode = 'never'
            new_schemes = self.mendel(self.initial_scheme.ids, self.initial_scheme.molecs,
                                            mutate_mols=self.params['start_changing_mols'], 
                                            n_population=self.params['init_population'],
                                            mutation_rate=self.params['initial_mutation_rate'],
                                            kill_rate=self.params['initial_kill_rate'], verbose=True,
                                            pop_str=pop_str, parent_key=self.params['initial_scheme_dir'].split('/')[1],
                                            smart_mode=smart_mode, curious_pool_ratio=self.params['curious_pool_ratio'],
                                            parent_loss=self.initial_scheme.loss, ensure_n_mols=True,
                                            forbid_crossover=True, boost_mols=False,
                                            max_curious_ids=self.params['max_curious_ids'])
        else:
            pop_str = self._pop_str
        
       
        
        
        for iteration in range(n_iterations):
            
            if not self.params['silent']: 
                start_time = datetime.now()
                print (self._run_counter, 'from', len(new_schemes), 'will run using PAC.', file=self._out_file)
                if self._hide_output: self._out_file.flush()
            self._run_pac(pop_str, silent=self.params['silent'], max_processes=max_processes)
            
            
            
            pop_str = self.params['pop_prefix']+str(iteration+1)
            self._pop_str=pop_str
            
            #Measuring performance
            pop_losses = []
            delta_losses = []
            #value = <value_if_true> if <expression> else <value_if_false>
            molecs_to_measure = self.imp_molecs if self.measure_only_imp else []
                
            for scheme_i in new_schemes:
                scheme_loss = scheme_i.get_loss(self._ideal_dfs, big_loss_def=self.params['big_loss_def'], weight_by=self.params['weight_by'], 
                                                imp_molecs=molecs_to_measure, mode=self.params['fitness_metric'], 
                                                 merge_fn=self.params['weight_fn'], #min_is_better=self.min_is_better, deprecated
                                                weight_factor=self.params['weight'], path=self.temp_path,
                                               abs_threshold = self.params['abs_threshold'])
                delta_loss = scheme_i.parent_loss - scheme_i.loss
                if scheme_loss!=None:
                    pop_losses.append(scheme_loss)
                    delta_losses.append(delta_loss)

                if (save_all_schemes and self._run_counter > 0) and not scheme_i._found:
                    self.master_df = scheme_i.save_scheme(self.master_df, templates=self._templates, 
                                                          imp_mols_flag=self.measure_only_imp, save_threshold=self.params['save_schemes_threshold'],
                                                         weights=self.params['weight'], photo_flag=self.params['photo_flag'],
                                                        single_planet_flag=self.params['single_planet_flag'], model_dict=self.params)
            if save_all_schemes:
                self.master_df.to_csv(self._master_df_path, index=False)
            if not self.params['silent']:
                print (len(pop_losses),'schemes survived.', file=self._out_file)
                if self._hide_output: self._out_file.flush()
                
            n_fittest = self.params['n_fittest']
            pop_losses = np.array(pop_losses)
            delta_losses = np.array(delta_losses)
            
            #  Reinitializing losses
            self._run_counter=0
            id_best_schemes = np.zeros((n_fittest), dtype=int)
            self.best_loss = self.params['big_loss_def']*np.ones((n_fittest))
            self.best_schemes = []

            if self.params['always_deterministic'] or iteration < self.params['n_before_wheel']:
                sorted_inds = np.argsort(pop_losses)
                if not self.min_is_better:
                    sorted_inds = sorted_inds[::-1]
                self.best_loss = pop_losses[sorted_inds[:n_fittest]]
                id_best_schemes = sorted_inds[:n_fittest]

                """
                for i in range(len(pop_losses)):
                    combined_loss  = pop_losses[i]
                    if np.any(combined_loss < self.best_loss):
                        self.best_loss = np.append(self.best_loss, [combined_loss])
                        id_best_schemes = np.append(id_best_schemes, [i])
                        inds = np.argsort(self.best_loss)
                        
                        self.best_loss = self.best_loss[inds]
                        id_best_schemes = id_best_schemes[inds]
                        self.best_loss = self.best_loss[:n_fittest]
                        id_best_schemes = id_best_schemes[:n_fittest]
                """

            else: #weighted roulette wheel
                sumfitness = np.sum(pop_losses)
                probabilities = [1 - ploss/sumfitness for ploss in pop_losses]
                inds_to_search = []
                for i in range(len(probabilities)):
                    prob = probabilities[i]
                    aux = i*np.ones(int(prob*len(probabilities)))
                    inds_to_search.extend(aux)

                for i_fit in range(n_fittest):
                    ind = np.random.randint(0, len(inds_to_search))
                    chosen_i = int(inds_to_search[ind])
                    id_best_schemes[i_fit] = chosen_i
                    self.best_loss[i_fit] = pop_losses[chosen_i]

            if self.params['add_per_delta'] > 0:
                n_add = self.params['add_per_delta']
                sorted_per_delta = np.argsort(delta_losses)[::-1]
                print (len(sorted_per_delta))
                k = 0
                _id = 0
                while (k < n_add and _id < len(sorted_per_delta)):
                    if not sorted_per_delta[_id] in id_best_schemes:
                        id_best_schemes = np.append(id_best_schemes, [sorted_per_delta[_id]])
                        k+=1   
                    _id+=1
                    

                    
                    
            for id_ in id_best_schemes:
                self.best_schemes.append(new_schemes[id_])
            
            if deallocate:    
                del new_schemes
                del pop_losses
                 
            if (np.any(np.max(self.best_loss) > self.old_losses) and iteration%self.params['mod_elitism']==0
                ) and iteration > 0:
                print('...Returning where you were happy...', file=self._out_file)
                if self._hide_output: self._out_file.flush()
                self.best_schemes.extend(self.old_schemes)
                self.best_loss=np.append(self.best_loss, self.old_losses)
                aux_best_schemes = []

                id_best_schemes = np.argsort(self.best_loss)[:n_fittest]
                self.best_loss = self.best_loss[id_best_schemes]

                for id_ in id_best_schemes:
                    sch = self.best_schemes[int(id_)]
                    aux_best_schemes.append(sch)
                self.best_schemes = aux_best_schemes
            
            self.best_loss_log.extend(self.best_loss)
            self.best_schemes_log.extend(self.best_schemes)
            
            self.optimate_losses=np.append(self.optimate_losses, [self.best_loss[0]])
            self.optimate_schemes.append(self.best_schemes[0])
                
            self.old_losses = self.best_loss
            self.old_schemes = self.best_schemes
            
            if iteration%n_print_loss_detail==0 and iteration > 0:
                if plot_loss:
                    plt.close('all')
                    fig = plt.figure(figsize=figsize)
                    ax = fig.add_subplot(111)
                    ax.plot(self.optimate_losses[1:], color='coral')
                    ax.set_ylabel(self.params['fitness_metric'])
                    ax.set_xlabel('Iterations')
                    if self.params['plt_show']:
                        plt.show()
                    fig.savefig('1D_PAC/loss_log.png', dpi=200)
                
                self.print_loss(self.best_schemes[0])
                self.best_schemes_dict = self.best_schemes_hist(self.best_schemes_log)
                

            if not self.params['silent']:
                print('.........................................................................................................', file=self._out_file)
                print ('Model losses: Min', np.min(self.best_loss), 'Max:', np.max(self.best_loss) , file=self._out_file)
                end_time = datetime.now()
                print ('Duration: {}'.format(end_time - start_time), file=self._out_file)
                print(file=self._out_file)
                if self._hide_output: self._out_file.flush()
            if not iteration + 1 == n_iterations or multiple_runs:
                n_population = self.params['n_population']
                crossover_round = iteration % self.params['n_crossover'] == 0

                if crossover_round and iteration > 1:
                    if not self.params['silent']:
                        print('                                       Crossover round, prepare for breeding.', file=self._out_file)
                        if self._hide_output: self._out_file.flush()
                    new_best_schemes, exp_max_comb = self._crossover(n_times=self.params['crossover_tries'], pop_str=pop_str,
                                                       max_comb=self.params['max_comb_crossover'])
                    if self._run_counter > 0:
                        print (self._run_counter, 'of ', exp_max_comb, 'will be calculated', file=self._out_file)
                        if self._hide_output: self._out_file.flush()
                        self._run_pac(pop_str=pop_str)
                    for init_sch in new_best_schemes:
                        _ = init_sch.get_loss(self._ideal_dfs, big_loss_def=self.params['big_loss_def'], 
                                                        weight_by=self.params['weight_by'], 
                                                imp_molecs=molecs_to_measure, mode=self.params['fitness_metric'], 
                                                merge_fn=self.params['weight_fn'], #min_is_better=self.min_is_better, deprecated
                                                weight_factor=self.params['weight'], path=self.temp_path)
                    n_population = crossover_pop
                else:
                    new_best_schemes = self.best_schemes
                    
                number_of_mols, number_of_reacs, imp_molecs_present, new_schemes = ([] for j in range(4))
 
                for initial_scheme in new_best_schemes:
                    imp_molecs_present.append(np.all([m in initial_scheme.molecs for m in self.imp_molecs]))
                    curious = iteration%self.params['mod_curious']==0
                    smart_mode = self.params['smart_choice_mode'] if not curious else 'curious'
                    
                    schemes = self.mendel(initial_scheme.ids, initial_scheme.molecs,
                                            mutate_mols=self.params['mutate_mols'], 
                                            n_population=n_population,
                                            mutation_rate=self.params['mutation_rate'],
                                            kill_rate=self.params['kill_rate'], verbose=False,
                                            pop_str=pop_str, start_id=len(new_schemes),
                                            parent_key=initial_scheme._fname_key,
                                            smart_mode=smart_mode,
                                            max_threshold_ids=self.params['max_smart_ids_per_reac'],
                                            smart_threshold=self.params['filter_ids_threshold'],
                                            parent_loss=initial_scheme.loss, ensure_n_mols=True,
                                            
                                            )

                    new_schemes.extend(schemes)
                    number_of_reacs.append(len(initial_scheme.ids))
                    number_of_mols.append(len(initial_scheme.molecs))
                     

                
                

                if not self.params['silent'] and len(number_of_mols) > 0:
                    print (file=self._out_file)
                    print ('New initial schemes:', file=self._out_file)
                    print ('number of molecs of them: ', np.min(number_of_mols),'-',np.max(number_of_mols), file=self._out_file)
                    print ('number of reacs of them:', np.min(number_of_reacs),'-',np.max(number_of_reacs), file=self._out_file)
                    print ('are all imp molecs there?', np.all(imp_molecs_present), file=self._out_file)
                    print (len(new_schemes), 'were created', file=self._out_file)
                    if self._hide_output: self._out_file.flush()
                   

            if ((iteration%n_write == 0 or iteration+1==n_iterations) or self.optimate_loss < self.params['stop_loss'] ) and iteration>min_iteration_write:
                
                
                self._write_schemes(reset_init=True, init_str=opt_str, write_all_best=not save_all_schemes, 
                                    loss_threshold=self.params['save_schemes_threshold'])
                self.print_loss(self.optimate_scheme)
                if save_all_schemes:
                    self.schemes_heatmap(apply_context=True, n_pops=iteration, 
                                     init_str=self.params['pop_prefix'], planets=self.params['planets'], 
                                     fname='1D_PAC/genes'+str(date.today())+'.pdf')
                
                if self.optimate_loss < self.params['stop_loss']:
                    print ('Stop condition fullfilled with loss:', self.optimate_loss, file=self._out_file)
                    if self._hide_output: self._out_file.flush()
                    break

if __name__ == "__main__":
    photochem = True
    small_run = False
    planets = ['hd2094', 'hd1897'] 

    pop_prefix = datetime.now().strftime("pop%H%M-")

    no_photo_dict = {
        'hd2094': ['1D_PAC/hd2094_extended.inp', '1D_PAC/hd2094_extended'],
        'hd1897': ['1D_PAC/hd1897_extended.inp', '1D_PAC/hd1897_extended'],
        'wasp107': ['1D_PAC/wasp107.inp', '1D_PAC/wasp107_full']
    }

    photo_dict = {
        'hd2094': ['1D_PAC/hd2094_photo.inp', '1D_PAC/hd2094_photo'],
        'hd1897': ['1D_PAC/hd1897_photo.inp', '1D_PAC/hd1897_photo'],
        'wasp107': ['1D_PAC/wasp107_photo.inp', '1D_PAC/wasp107_photo']
    }

    templates = []
    full_dirs = []
    for pl in planets:
        if photochem:
            templates.append(photo_dict[pl][0])
            full_dirs.append(photo_dict[pl][1])
        else:
            templates.append(no_photo_dict[pl][0])
            full_dirs.append(no_photo_dict[pl][1])

    pop_prefix = datetime.now().strftime("pop%H%M-")

    if photochem:
        initial_dir = '1D_PAC/opt0705_photo'   
    else:
        initial_dir = '1D_PAC/pca_alone46'

    if small_run:
        init_population, n_fittest, n_population, n_iterations = 25, 5, 5, 5
    else:
        init_population, n_fittest, n_add, n_population, n_iterations = 500, 25, 2, 15, -1

    GA = Darwen(mutation_rate=0.1, mutation_prob =.2, crossover_prob= 0.4, kill_rate=0.15, 
                mod_elitism=7, n_crossover=1e3, 
                initial_mutation_rate=0.2, initial_kill_rate=0.1, 
                molecs_low_limit=30, molecs_high_limit=50, weight=[1e-2, 1],  
                pop_prefix=pop_prefix, search_schemes_first=True, planets=planets,
                init_population=init_population, n_fittest=n_fittest, n_population=n_population,
                initial_scheme_dir=initial_dir, full_dirs=full_dirs, templates=templates, 
                max_curious_ids=3, update_ids_df=True, rosetta_path='1D_PAC/opt',
                hide_output=True, photo_flag=photochem, single_planet_flag=len(planets)==1,
                weight_by='molecs&logmajor', master_df_name='master_index.csv',
                extended_range_flag=False, add_per_delta=n_add)
    GA.clean_dir()
    GA.run(n_iterations=-1, n_write=3, save_all_schemes=True, max_processes=40)