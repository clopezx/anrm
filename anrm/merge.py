#merge.py

from pysb import *

class Edit_Monomers(object):
    def __init__(self, all_monomers, edits):
        self.less_monomers   = self.del_edit_monomers(all_monomers, edits)
        self.merged_monomers = self.merge_monomers(self.less_monomers, edits)

    def check_duplicate_monomers(self, all_monomers):
        list_models = all_monomers.keys()
        duplicates = []
        for i in range(len(list_models)):
            for j in all_monomers[list_models[i]]:
                dups_models = []
                if i+1 < len(list_models):
                    for k in range(len(list_models))[i+1:]:
                        if j.name in all_monomers[list_models[k]]._map:
                            dups_models.append(list_models[k])
                    if len(dups_models) > 0:
                        duplicates.append(j.name)
                        dups_models.append(list_models[i])
                        print "Found duplicate monomer: ", j.name,
                        print " in the following models: ", dups_models
        
        return duplicates

    def del_edit_monomers(self, all_monomers, edits):
        #remove duplicate monomers from all of the models' ComponentSets
        less_monomers = all_monomers
        monomer_names = all_monomers.keys()
        for j in edits.keys():
            monomer_index = less_monomers._index_map[j]
            to_remove = less_monomers[monomer_index:monomer_index+1]
            less_monomers = less_monomers - to_remove
       
        return less_monomers
    
    def merge_monomers(self, less_monomers, new):
        from pysb.core import ComponentSet
        # I should include a verification step that all duplicates are accounted for. Are all duplicates accounted for? are there any extra monomers?
        def create_monomers(new):
            Model('merger')
            for name in new.keys():
                Monomer(name, new[name][0], new[name][1])
            return merger.monomers
    
        merged_monomers = ComponentSet()
        new_monomers = create_monomers(new)
        
        merged_monomers |= less_monomers
        merged_monomers |= new_monomers
    
        return merged_monomers

class Edit_Rules(object):
    def __init__(self, all_rules, edits):
        self.less_rules   = self.del_edit_rules(all_rules, edits)
        self.merged_rules = self.merge_rules(self.less_rules, edits)

    def del_edit_rules(self, all_rules, edits):
        #remove duplicate monomers from all of the models' ComponentSets
        less_rules = all_rules
        rule_names = all_rules.keys()
        for j in edits.keys():
            rule_index = less_rules._index_map[j]
            to_remove = less_monomers[rule_index:rule_index+1]
            less_rules = less_rules - to_remove
        
        return less_rules

    def merge_rules(self, less_monomers, new):
        from pysb.core import ComponentSet
        # I should include a verification step that all duplicates are accounted for. Are all duplicates accounted for? are there any extra monomers?
        def create_rules(new):
            Model('merger_rules')
            from irvin_AN_crosstalk import rule_revisions
            rule_revisions
            return merger.rules
        
        merged_rules = ComponentSet()
        new_rules = create_rules(new)
        
        merged_rules |= less_rules
        merged_rules |= new_rules
        
        return merged_rules