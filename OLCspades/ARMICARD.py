import time
from collections import defaultdict

__author__ = 'mike knowles'


class Card:
    """
    CARD requires a gene # as an input class has three functions:
    Card(antidict, gene)
    resist will :return and list of antibiotics for a relevant gene
        Includes functionality to trace the dependencies of a gene complex
        Utilizes recursion to achieve all possible antibiotics
        .resist(genome)
    anti will :return a list of antibiotics for a given gene
        .anti()
    sens will :return a list of sensitivities for a given gene
        .sens()
    """

    def __init__(self, antidict, index, plusdict=None):  # Initialize class and inherit self
        self.index = index  # Defines a gene number that can be passed to functions within the class
        self.antidict = antidict
        self.plusdict = plusdict

    def resist(self, genome=None, gene=None, tolc=None):  # Begin resist function and import initialized self
        resistlist = []  # Initialize dict
        genedict = self.antidict[self.index]
        deps = True
        keystr = self.index.encode('utf8') + ", " + genedict['name'].encode('utf8')
        if "member" in genedict and genome is not None:  # check if dependencies are satisfied
            deps = False
            membstr = lambda x: x.encode('utf8') + ", " + self.antidict[x]['name'].encode('utf8')
            index = {keystr: [membstr(memb) for memb in genedict['member']
                              if memb in self.plusdict and memb != tolc]}
            if len(index[keystr]) == len(genedict['member']):
                # check if the list of requirements is complete
                deps = True
                index[keystr].sort()
                # resistlist.extend([dict((resist, mdict) for resist in genedict['resist'])])  # create list of dicts
            else:
                index = []
        else:
            index = {keystr: gene} if gene else [keystr]
        if "resist" in genedict and deps:  # If the key "resist" in gene
            if "complex" in genedict and genome is not None:
                '''If the key complex in gene defines their are depenedenies'''
                # count = 0  # for each resistance set count at zero
                for comp in genedict["complex"]:  # Allow for multiple dependencies
                    inde = Card(self.antidict, comp, self.plusdict).resist(genome, tolc=tolc)
                    if inde:
                        resistlist.extend(inde)
                    # recurse through the same class if complexes are satisfied extend the list
            else:  # if no complex then just return the list
                resistlist.extend([dict((resist, index) for resist in genedict['resist'])])
        if "isa" in genedict and deps:  # Recursion for parent antibiotic traits
            for depend in genedict["isa"]:
                for amr in Card(self.antidict, depend, self.plusdict).resist(genome, index, tolc=tolc):
                    # Call self to recurse through the same class
                    # if amr not in resistlist:
                    resistlist.append(amr)
                    # resistlist.extend(self.resist(genome))  # Call self to recurse through the same class
        return resistlist  # return the extended list

    def anti(self):
        if "resist" in self.antidict[self.index]:
            return self.antidict[self.index]["resist"]

    def sens(self):
        if "sensitivity" in self.antidict[self.index]:
            return self.antidict[self.index]["sensitivity"]

    def function(self):
        if "function" in self.antidict[self.index]:
            return self.antidict[self.index]["function"]


class DictBuild:
    """Simple class to build a list or dictionary without repeats"""

    def __init__(self, index):
        self.key = index

    def add(self, lst):
        if not self.key:
            self.key = lst
        else:
            for drug in lst:
                if drug not in self.key:
                    self.key.append(drug)
        # self.key = sorted(self.key, key=lambda s: s.lower())
        return self.key


def recur(current, existing, index, dup=True):
    for item in current:
        if item in existing:
            citem, eitem = current[item], existing[item]
            if type(citem) is list:
                if type(eitem) is list and citem[0] not in eitem:
                    if len(citem) > 1:
                        clist = current
                    else:
                        clist = citem[0]
                    existing[item].append(clist)
                elif citem[0] not in eitem:
                    existing[item][item] = citem
                elif citem[0] in eitem:
                    return False, existing
            elif type(eitem) is list:
                eitem.append(citem)
            else:
                dup, existing[item] = recur(citem, eitem, index+1, dup)
        elif not type(existing) in (str, unicode, list) and index != 0:
            existing[item] = current[item]
    return dup, existing


def decipher(plusdict, antidict, outputs, metadata, tolc=None):
    from copy import deepcopy
    import json
    outputdict = {}
    for genome in sorted(plusdict):  # iterate through plus dict
        outputdict[genome] = {"resist": defaultdict(list), "sensitivity": [], "genes": []}
        arodi = defaultdict(list)
        resistance = outputdict[genome]["resist"]
        for gene in plusdict[genome]:
            if 'name' in antidict[gene]:
                plusdict[genome][gene].insert(0, antidict[gene]['name'])
            analysis = Card(antidict, gene, plusdict[genome])
            if plusdict[genome][gene]:
                sens = analysis.sens()  # check sensitivities
                for resist in analysis.resist(genome, tolc=tolc):  # check resistances
                    if resist is not None:
                        for aro in resist:
                            if resist[aro] not in arodi[aro]:
                                arodi[aro].append(resist[aro])
                                if type(resist[aro]) is dict:
                                    cpa = deepcopy(resistance[aro])
                                    new, dup = [], True
                                    for existing in cpa:
                                        up, newitem = recur(resist[aro], existing, 0)
                                        new.append(newitem)
                                        dup = up if not up else dup
                                    if new == resistance[aro] and dup:
                                        resistance[aro].append(resist[aro])
                                    elif new and dup:
                                        resistance[aro] = new
                                elif resist[aro][0] not in resistance[aro]:
                                    resistance[aro].extend(resist[aro])
                                resistance[aro].sort()

                outputdict[genome]["genes"].append(tuple([gene] + plusdict[genome][gene]))
                if sens is not None:
                    outputdict[genome]["sensitivity"] = DictBuild(outputdict[genome]["sensitivity"]).add(sens)
        outputdict[genome]["genes"].sort(key=lambda tup: tup[0])
        outputdict[genome]["sensitivity"].sort()
    json.dump(outputdict,
              open("%s/ARMI_CARD_results.json" % outputs, 'w'),
              sort_keys=True,
              indent=4,
              separators=(',', ': '))
    antilist = []

    for gene in antidict:  # build hearder list
        resistances = Card(antidict, gene).anti()
        if resistances is not None:
            for resist in resistances:
                if resist not in antilist:
                    antilist.append(resist)
    antihead = "Genome"
    drugcounter = {}
    antilist = sorted(antilist, key=lambda s: s.lower())  # sort header case insensitive
    for anti in antilist:
        antihead += ",\"%s\"" % anti
        drugcounter[anti] = 0

    antistr = ""

    ''' Build csv string '''

    # for genome in sorted(outputdict):

    for sample in metadata:
        # genomename = '\n{}'.format(os.path.split(os.path.splitext(genome)[0])[1].replace('_filteredAssembled', ""))
        # print genome, genomename
        # genomename = '\n' + sample.name
        antistr += '\n' + sample.name
        genomecount = 0
        for drug in antilist:
            if drug in outputdict[sample.name]["resist"]:
                antistr += ",%i" % len(outputdict[sample.name]["resist"][drug])
                drugcounter[drug] += 1
                genomecount += 1
            else:
                antistr += ",-"
        antistr += ",%i" % genomecount
        # Open the report
        if sample.general.bestassemblyfile != 'NA':
            with open('{}{}_{}.csv'.format(sample['ARMI'].reportdir, sample.name, 'ARMI'), 'wb') as report:
                # Write the row to the report
                report.write(antihead)
                report.write(antistr)

    antihead += "\nCount"
    for drug in antilist:
        antihead += ",%i" % drugcounter[drug]
    antihead += antistr
    with open("%s/ARMI_CARD_results.csv" % outputs, 'w') as f:
        f.write(antihead)
