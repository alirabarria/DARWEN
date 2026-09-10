from darwen_functions import *
import warnings

from datetime import datetime, date
class Scheme:
    """
    Represents a chemical reaction scheme, encapsulating the functionality to search, save, and evaluate schemes
    based on their performance against ideal data. This class is integral to optimizing chemical reaction paths
    in a given chemical space, leveraging a genetic algorithm.
    """
    def __init__(self, ids, molecs, master_df, parent_key='', scheme_id=0, out_path='1D_PAC/gen/', pop_str='pop0', 
                 planets=['hd2094', 'hd1897'], parent_loss=1e5, len_reacs=958, max_mols=107,
                 len_ids_pool=-1, schemes_path='saved_schemes/', search_first=True, change_fname_key=True,
                 curiosity=0, reacs_dict={}, specs_from_file={}, specs_lines=[]):
        """
        Initializes a Scheme object with given identifiers, molecules, and configuration parameters.

        Parameters:
        - ids (array-like): Identifiers for the reactions involved in the scheme.
        - molecs (array-like): Molecules involved in the scheme.
        - master_df (DataFrame): The master dataframe containing all schemes' information.
        - parent_key (str): Identifier for the parent scheme, if any.
        - scheme_id (int): Unique identifier for this scheme.
        - out_path (str): Output directory for saving scheme files.
        - pop_str (str): Population string identifier for scheme categorization.
        - planets (list of str): Planets for which the scheme is evaluated.
        - len_ids_pool (int): The length of the identifiers pool.
        - schemes_path (str): Directory where schemes are saved.
        - search_first (bool): Whether to search for the scheme in the existing database initially.
        """
        self.ids = np.sort(np.unique(ids))
        self._len_genes=len_reacs
        self.genes = self.ids_to_genes(self.ids)
        self.molecs = np.sort(np.unique(molecs))
        self.loss, self._fname_key, self._found, self._found_data, self._loss_type = None, '', False, [], ''
        self._scheme_id, self._pop_str, self._out_path, self._planets = scheme_id, pop_str, out_path, planets
        self._len_ids_pool, self._schemes_path, self._search_first = len_ids_pool, schemes_path, search_first
        self._parent_key, self.parent_loss, self._curiosity, self._len_reacs =parent_key, parent_loss, curiosity, len_reacs
        self._reacs_dict = reacs_dict
        self._specs_from_file = specs_from_file
        self._specs_lines = specs_lines

        
        if search_first:
            self._found, self._found_path = self.search_scheme(self._schemes_path, master_df, change_fname_key=change_fname_key)

    def ids_to_genes(self, reacs_array):
        chromosome = np.zeros((self._len_genes), dtype=bool)
        inds = np.array(reacs_array).astype(int)
        chromosome[inds] = True
        return chromosome
    
    def genes_to_ids(self):
        self.ids = np.where(self.genes)[0]

    def self_mutate(self, mutation_prob=.01):
        len_genes = self._len_genes
        ids_to_change = np.random.randint(0, len_genes, size=int(len_genes*mutation_prob))
        genes_to_change = self.genes[ids_to_change]
        self.genes[ids_to_change] = [not aux for aux in genes_to_change]
        self.genes_to_ids()

    def search_scheme(self, directory, master_df, change_fname_key=True):
        """
        Searches for an existing scheme within a specified directory, comparing against a master DataFrame that contains metadata for all existing schemes.

        Parameters:
        - directory (Path): The directory to search in, specified as a Pathlib object for OS-independent path resolution.
        - master_df (DataFrame): A pandas DataFrame containing metadata for existing schemes. Expected to have columns for 'len_pool', 'mols_len', 'ids', and 'directory' at minimum, which are used to identify matching schemes based on molecule counts and identifier ranges.

        Returns:
        - found (bool): True if a matching scheme is found, False otherwise.
        - path_to_scheme (Path): The Path object pointing to the found scheme's directory. If not found, this will point to the intended directory based on the scheme's specifications.
        """
        directory = Path(directory) if not isinstance(directory, Path) else directory
        len_mols = len(self.molecs)
        directory = directory / Path(str(int(len_mols)))
        ids = self.ids
        min_id, mean_id, max_ids = int(np.min(ids)), int(np.mean(ids)), int(np.max(ids))
        scheme_directory = Path(str(int(len_mols))+'_'+str(min_id)+'-'+str(mean_id)+'-'+str(max_ids))
        path_to_scheme = directory / scheme_directory

        if len(master_df)>0:
            mols_df = master_df[master_df.mols_len == len(self.molecs)]
            mols_df = mols_df[mols_df.ids_len == len(self.ids)]
            if self._len_ids_pool > -1:
                existing_df = mols_df[mols_df.len_pool == self._len_ids_pool]
            else:
                existing_df = mols_df
                    
            ids_condition = [np.all(from_df_to_array(aux)== ids)  for aux in existing_df.ids]
            ids_df = existing_df[ids_condition]
            if len(ids_df) > 1:
                warnings.warn("Not unicity of ids.", UserWarning)
            if len(ids_df)>0:
                self._found = True
                row = ids_df.iloc[0]
                if row['directory'] != 'Not saved' and change_fname_key:
                    self._fname_key = row['fname_key']
                    self._found_data = row.to_dict()
 
            if not self._found:
                mask_dirs = [str(path_to_scheme) in aux for aux in mols_df['directory']]
                direcs_df = mols_df[mask_dirs]
                #mask_dirs = [path_to_scheme in _dir for _dir in master_df['directory']]
                if len(direcs_df) > 0:
                    aux=str(path_to_scheme)+'_'+str(len(direcs_df))
                    path_to_scheme = Path(aux)
                #else path_to_scheme = itself
 
        return self._found, path_to_scheme
        
    
    def save_scheme(self, master_df, out_path='', templates=[], imp_mols_flag=False, save_threshold=3,
                   weights=[], single_planet_flag=False, photo_flag=False, merge_fn=np.max, force_save=False,
                   change_fname_key=True, directories=[], model_dict={}, force_ind=-1):
        """
        Saves the scheme to disk if it meets certain criteria, such as having a loss below a specific threshold, and updates the master DataFrame with the scheme's details.

        Parameters:
        - master_df (DataFrame): The master DataFrame to update with the scheme's details upon saving.
        - out_path (Path): The base path where the scheme should be saved, specified as a Pathlib object.
        - templates (list): A list of template files or strings used in generating the scheme files.
        - imp_mols_flag (bool): A flag indicating whether to consider only important molecules in the scheme's evaluation.
        - save_threshold (float): The loss threshold below which the scheme is considered worth saving.
        - weights (list): A list of weights to apply to different aspects of the scheme during evaluation.

        Returns:
        - DataFrame: The updated master DataFrame with the scheme's details if saved.
        """
        
        if out_path=='':
            out_path=self._out_path
            

        found, path_to_scheme = self.search_scheme(self._schemes_path, master_df, change_fname_key=change_fname_key)
        if len(directories) > 0:
            mask = [str(path_to_scheme) in aux for aux in directories]
            repeated = np.sum(mask)
            if repeated > 0:
                path_to_scheme = Path(str(path_to_scheme) + '_' + str(int(repeated))).resolve()

        if not found or force_save:
            if merge_fn(self.losses) < save_threshold or force_save:
                path_to_scheme.parent.mkdir(parents=True, exist_ok=True)
                move_files_with_string(str(self._out_path), str(path_to_scheme), self._fname_key+'_', copy=True)
                self._write_scheme(out_path=path_to_scheme, write_ids=True, templates=templates, imp_mols_flag=imp_mols_flag, check_dir=False)
            else:
                path_to_scheme = 'Not saved'
            delta_loss = self.parent_loss - self.loss
            if force_ind == -1:
                ind_df = len(master_df)
            else:
                ind_df = force_ind
            self.genes_list = self.genes.astype(int).tolist()
            #                      ['mols_len', 'ids_len' 'ids', 'molecs',   
            master_df.loc[ind_df]=[len(self.molecs), len(self.ids), self.ids.tolist() , self.molecs.tolist(), 
                                    # 'loss', 'delta_loss', 'losses',       
                                    self.loss, delta_loss, self.losses, 
                                    # 'abs_losses', 'all_rel_losses', 'major_loss',  'directory',     
                                    self.max_err, self._rel_all, self._major_loss, str(path_to_scheme),
                                    #'fname_key', 'parent', 'len_pool', 
                                    self._fname_key,  self._parent_key, self._len_ids_pool, 
                                    # 'date', 'loss_type', 'worst_mols', 'weights', 'planets'
                                    str(date.today()), self._loss_type, self._worst_mols, weights, self._planets,
                                    # 'photo_flag', 'single_planet_flag', 'curiosity', 'genes',
                                    photo_flag, single_planet_flag, self._curiosity, self.genes_list,
                                    # 'model_dict']
                                    model_dict]
        if len(directories) > 0:
            directories.append(str(path_to_scheme))
            return master_df, directories
        else:
            return master_df
            
    
        
    def _write_scheme(self, write_ids=False, check_dir=True, init_str='', imp_mols_flag=False, out_path='', templates=[], 
                      change_spec_file=True):
        """
        A helper method to write scheme details to files, including creating necessary input files for simulations based on the scheme.

        Parameters:
        - write_ids (bool): If True, includes the scheme's identifiers in the output files.
        - check_dir (bool): If True, checks if the output directory exists and creates it if not.
        - init_str (str): An initial string or header to include at the top of the output files.
        - imp_mols_flag (bool): If True, focuses on important molecules in the output files.
        - out_path (Path): The base path for output files, specified as a Pathlib object.
        - templates (list): A list of template files or strings used as bases for the output files.

        Returns:
        - None: The method writes files to disk based on the scheme and does not return a value.
        """
        if out_path=='':
            out_path = Path(self._out_path) if not isinstance(self._out_path, Path) else self._out_path
        if init_str=='':
            init_str=self._pop_str

        reacs_dict = self._reacs_dict
        if check_dir:
            (out_path / Path(init_str)).mkdir(parents=True, exist_ok=True)

        if not self._search_first or not self._found:
            self._fname_key=init_str+'_'+str(self._scheme_id)

        for id_planet in range(len(self._planets)):
            planet_flag = self._planets[id_planet]
            fname = self._fname_key
            
            inp_lines= templates[id_planet]
            
            with open(out_path / Path(fname+'_'+planet_flag+'.inp'), 'w') as new_inp_file:
                if change_spec_file and len(self._specs_lines) > 0:
                    inp_lines[0] = fname+'.spec              ! species file\n'
                inp_lines[3] = fname+'.reac               ! reaction file\n'
                for line in inp_lines:
                    new_inp_file.write(line)
               
            pac_fname = init_str+'/point_to_'+fname+'_'+planet_flag+'.inp'
            if check_dir:
                with open(out_path / Path(pac_fname), 'w') as new_pac_file:
                    new_pac_file.write(fname+'_'+planet_flag+'.inp')
            
            if not self.loss == None:
                write_new_scheme(self.ids, self.molecs, fname=str(out_path)+'/'+fname, write_ids=write_ids, loss=self.losses[id_planet],
                                imp_mols_flag=imp_mols_flag, reacs_dict=reacs_dict, specs_lines=self._specs_lines,
                                specs_from_file=self._specs_from_file)
            else:
                write_new_scheme(self.ids, self.molecs, fname=str(out_path)+'/'+fname, write_ids=write_ids,
                                imp_mols_flag=imp_mols_flag, reacs_dict=reacs_dict, specs_lines=self._specs_lines,
                                specs_from_file=self._specs_from_file)


                
    def get_loss(self, ideal_dfs_array, big_loss_def=1e5, weight_by='molecs', imp_molecs=[], 
                 mode='max_percent', merge_fn=np.mean, weight_factor=[1e-3, 1e6], path='',
                abs_threshold=1e5, major_threshold=1e-20, curiosity=0, force_run=False,
                error_verbose=False, ideal_dirs=[], insert_dirs=[]):
        """
        Calculates the loss of the scheme by comparing it against ideal datasets. The loss calculation can be adjusted by various parameters to focus on specific molecules, use different metrics, and apply weighting.

        Parameters:
        - ideal_dfs_array (list of DataFrame): A list of pandas DataFrames representing ideal datasets for comparison.
        - big_loss_def (float): A default large loss value used as a penalty in certain conditions.
        - weight_by (str): Determines how the loss is weighted; options include 'molecs' for molecule-based weighting.
        - imp_molecs (list of str): Specifies important molecules to focus on in the loss calculation.
        - mode (str): The metric used for calculating loss ('max_percent', 'mre', etc.).
        - merge_fn (function): A function used to aggregate multiple loss values, typically np.mean or np.max.
        - weight_factor (list of float): Factors used to weight different aspects of the loss calculation.
        - path (Path): The path to additional resources needed for loss calculation, specified as a Pathlib object.
        - abs_threshold (float): An absolute threshold for loss beyond which the scheme is considered unviable.
        - calculate_for_all (bool): If True, calculates loss for all provided ideal datasets.

        Returns:
        - float: The calculated loss for the scheme.
        """
        
        self.losses=[]
        self.max_err=[]
        self._worst_mols=[]
        self._rel_all = []
        self._major_loss = []
        self._dfs_to_return = []
        self._loss_type=mode+'-'+weight_by
        if curiosity != 0:
            self._curiosity = curiosity
        if path=='':
            path = self._out_path
            
        all_ready=False
        
        if self._found and not force_run:
            path = self._found_path
            row = self._found_data
            if row['loss_type']==self._loss_type and np.all(from_df_to_array(row['planets'])==self._planets):
                all_ready=True
                self.loss=row['loss']
                self.losses=from_df_to_array(row['losses'])
                self.max_err=from_df_to_array(row['abs_losses'])
                self._worst_mols=from_df_to_array(row['worst_mols'])

                self._rel_all=from_df_to_array(row['all_rel_losses'])
        if len(ideal_dirs) > 0:
            ideal_dfs_array = []
            for ideal_dir in ideal_dirs:
                ideal_df, _, _ = out_to_dfs(dir_=ideal_dir, error_purposes=True, return_only_mols=False)
                ideal_dfs_array.append(ideal_df)


        if not all_ready or len(insert_dirs) > 0:
            if len(ideal_dirs) + len(insert_dirs) > 0:
                aux_dfs = []

            for ideal_i in range(len(ideal_dfs_array)):
                ideal_df = ideal_dfs_array[ideal_i]
                if len(ideal_dirs) > 0:
                    aux_dfs.append(ideal_df)

                if len(insert_dirs) > 0:
                    inserted_dir = insert_dirs[ideal_i]
                    df_to_test, cons_err, _err = out_to_dfs(dir_= inserted_dir, error_purposes=True, return_only_mols=False)
                    aux_dfs.append(df_to_test)
                else:
                    df_to_test, cons_err, _err = out_to_dfs(dir_=path, error_purposes=True, file_key=self._fname_key+'_'+self._planets[ideal_i]+'_out', 
                                                  return_only_mols=True)
                    

                err_pop, max_err_pop, mol_pop, max_layer, mols_errors = calculate_error(df_to_test, ideal_df, imp_molecs=imp_molecs, punish_cons_err=cons_err,
                                                            mode=mode, big_loss_def=big_loss_def, major_threshold=major_threshold,
                                                            verbose=error_verbose)
                self.losses.append(err_pop)
                self.max_err.append(max_err_pop)
                self._worst_mols.append([mol_pop, max_layer])

                if len(insert_dirs) + len(ideal_dirs) > 0:
                    self._dfs_to_return.append(aux_dfs)

                if not len(imp_molecs)==0 and mols_errors != {}:
                    max_err_all = 0
                    major_loss = 0
                    for m in mols_errors.keys():
                        if mode == 'max_percent':
                            max_err_all = np.max([mols_errors[m][0], max_err_all])
                            major_loss = np.max([mols_errors[m][2], major_loss])

                        elif mode == 'mre':
                            max_err_all+=mols_errors[m][0]
                            major_loss+=mols_errors[m][2]
                    if mode=='mre':
                        max_err_all = max_err_all/len(mols_errors)
                else:
                    max_err_all = err_pop
                    major_loss = err_pop
                self._rel_all.append(max_err_all)
                self._major_loss.append(major_loss)

            if weight_by=='molecs':
                to_add = weight_factor[0]*len(self.molecs)
            elif weight_by=='molecs&allabs':
                abs_weighted = weight_factor[1] * merge_fn(self.max_err)
                if abs_threshold < abs_weighted:
                    to_add = weight_factor[0] * len(self.molecs)
                else:
                    to_add = weight_factor[0] * len(self.molecs) + abs_weighted
            elif weight_by=='molecs&relmajor':
                rel_weighted = weight_factor[1]*merge_fn(self._major_loss)
                to_add = weight_factor[0] * len(self.molecs) + rel_weighted
            elif weight_by=='molecs&logmajor':
                rel_weighted = weight_factor[1]*np.log10(merge_fn(self._major_loss))
                to_add = weight_factor[0] * len(self.molecs) + rel_weighted
            elif weight_by=='molecs&logall':
                all_weighted = weight_factor[2]*np.log10(merge_fn(self._rel_all))  
                rel_weighted = weight_factor[1]*np.log10(merge_fn(self._major_loss))
                to_add = weight_factor[0] * len(self.molecs) + rel_weighted + all_weighted
            self.loss=merge_fn(self.losses)+to_add-curiosity
        if len(insert_dirs)==0:
            return self.loss
        else:
            return self.loss, self._dfs_to_return