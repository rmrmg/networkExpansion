import json, argparse
from sys import exit
from itertools import product
from os.path import isfile
from rdkit import Chem
from rdkit.Chem import AllChem


class WrongSmarts(Exception):
    pass


class WrongSmiles(Exception):
    pass


class DuplicatedInputData(Exception):
    pass


def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--substrates', type=str, required=True, help='file with substrates in SMILES format')
    parser.add_argument('-g', '--generations', type=int, default=1, help='number of synthetic generation')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose mode - additional information will be print to stdout')
    parser.add_argument('-o', '--output', type=str,
                        help='name of output file if not given results will be printed to stdout')
    parser.add_argument('-r', '--reactions', type=str, required=True,
                        help='file with reactions database, see comment in code about expected format')
    args = parser.parse_args()
    return args


class Reactor:
    def __init__(self):
        self.mols = []
        self.smiles = set()
        self.reactionMatches = dict()
        self.reactionDoneCombinations = dict()

    @staticmethod
    def _getMatchesToRx(mol, rxn):
        for incogr in rxn.incompatGroups:
            if mol.HasSubstructMatch(incogr):
                return []
        pozs = []
        for poz, core in enumerate(rxn.reactants):
            if mol.HasSubstructMatch(core):
                pozs.append(poz)
        return pozs

    @staticmethod
    def _getMatchingCombination(allMatchesListDict):
        itr = product(*[allMatchesListDict[x] for x in allMatchesListDict])
        return itr

    def _addMatchesForRx(self, rxn, mols, begIdx):
        for idx, mol in enumerate(mols):
            matches = self._getMatchesToRx(mol, rxn)
            if not matches:
                continue
            if rxn.idx not in self.reactionMatches:
                self.reactionMatches[rxn.idx] = {idx: [] for idx in range(len(rxn.reactants))}
            for poz in matches:
                self.reactionMatches[rxn.idx][poz].append(begIdx + idx)

    def addCompounds(self, mols, rxes):
        begIdx = len(self.mols)
        for mol in mols:
            smiles = Chem.MolToSmiles(mol)
            if smiles in self.smiles:
                print("duplicated compounds")
                raise WrongSmiles
            self.smiles.add(smiles)
            self.mols.append(mol)
        for rxn in rxes:
            self._addMatchesForRx(rxn, mols, begIdx)

    def makeGeneration(self, rxes):
        products = []
        for rx in rxes:
            if rx.idx not in self.reactionMatches:
                continue
            for molIdxList in self._getMatchingCombination(self.reactionMatches[rx.idx]):
                if rx.idx in self.reactionDoneCombinations:
                    if tuple(molIdxList) in self.reactionDoneCombinations[rx.idx]:
                        continue
                else:
                    self.reactionDoneCombinations[rx.idx] = set()
                self.reactionDoneCombinations[rx.idx].add(tuple(molIdxList))
                molecules = [self.mols[idx] for idx in molIdxList]
                products.extend(rx.performRx(molecules))
        return products

    def filterOutSeenCompounds(self, molList):
        mols = []
        addedSmiles = set()
        for mol in molList:
            smiles = Chem.MolToSmiles(mol)
            if smiles in self.smiles or smiles in addedSmiles:
                continue
            addedSmiles.add(smiles)
            mols.append(mol)
        return mols


class Reaction:
    def __init__(self, jsontext):
        data = json.loads(jsontext)
        self.rxsmarts = data['reactionSmarts']
        self.rxn = AllChem.ReactionFromSmarts(data['reactionSmarts'])
        self.idx = data['idx']
        if not self.rxn:
            print(f"incorrect reaction smarts: {data['reactionSmarts']}")
            print("please fix your reaction database")
            raise WrongSmarts
        self.incompatGroups = self._getMolListFromSmartsList(data['incompatGroups'])
        self.bannedGroups = self._getMolListFromSmartsList(data["bannedProducts"])
        self.reactants = tuple(self.rxn.GetReactants())

    @staticmethod
    def _getMolListFromSmartsList(incolist):
        molList = []
        smaSet = set()
        for sma in incolist:
            if sma in smaSet:
                print(f"smarts {sma} is duplicated")
                print("please fix your reaction database")
                raise DuplicatedInputData
            mol = Chem.MolFromSmarts(sma)
            if not mol:
                print(f"incorrect smarts {sma}")
                print("please fix your reaction database")
                raise WrongSmarts
            molList.append(mol)
        return molList

    def getMatches(self, mol):
        matchedPoz = []
        for patt in self.incompatGroups:
            if mol.HasSubstructMatch(patt):
                return []
        for poz, patt in enumerate(self.reactants):
            if mol.HasSubstructMatch(patt):
                matchedPoz.append(poz)
        return tuple(matchedPoz)

    def performRx(self, mols):
        prods = self.rxn.RunReactants(mols)
        prods = [p[0] for p in prods]
        try:
            _ = [Chem.SanitizeMol(p) for p in prods]
        except:
            print("cannot sanitize one of reaction product probably due to incorrect reaction smarts")
            raise WrongSmarts
        return prods


def loadSubstrates(fn):
    # load substrates from file and return as Rdkit mol objects
    fh = open(fn)
    mols = []
    canonSmiles = set()
    for line in fh:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        for smiles in line.split('.'):
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                raise WrongSmiles
            canon = Chem.MolToSmiles(mol)
            if canon in canonSmiles:
                print(f'molecule {smiles} which has canonical form {canon} is duplicated in substrates')
                print("please correct your input data")
                raise DuplicatedInputData
            canonSmiles.add(canon)
            mols.append(mol)
    fh.close()
    return mols


def loadReactions(fn):
    # load reaction from text where each line contains information about reaction in json format
    # each line/json need to have following keys:
    # - idx: unique identificator of reaction, can be of any hashable type
    # - rxSmarts: chemical reaction in SMARTS format,
    # - incompatGroups: functional groups which are not allow in substrates
    # - bannedProducts: structural motifs which are not allow in product
    # Please note: quality of results *strongly* depend on proper definition
    #  of reaction SMARTS or SMARTS of functional groups.
    reactions = []
    fh = open(fn)
    for line in fh:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        rx = Reaction(line)
        reactions.append(rx)
    fh.close()
    return reactions


def printResultToFile(molList, genNum, fout):
    for mol in molList:
        smiles = Chem.MolToSmiles(mol)
        if fout is None:
            print(genNum, smiles)
        else:
            print(genNum, smiles, file=fout)


if __name__ == "__main__":
    args = parseArgs()
    substrates = loadSubstrates(args.substrates)
    if args.verbose:
        print(f"loaded {len(substrates)} substrates")
    reactions = loadReactions(args.reactions)
    if args.verbose:
        print(f"loaded {len(reactions)} reactions")
    reactor = Reactor()
    reactor.addCompounds(substrates, reactions)
    results = dict()
    if args.output and isfile(args.output):
        print("file exists remove/rename it or use other name")
        exit(1)
    fout = None
    if args.output:
        fout = open(args.output, 'w')
    for gennum in range(1, args.generations + 1):
        products = reactor.makeGeneration(reactions)
        if args.verbose:
            print(f"total products in {gennum} generation: {len(products)}")
        newProducts = reactor.filterOutSeenCompounds(products)
        if args.verbose:
            print(f"total unique products in {gennum} generation: {len(newProducts)}")
        printResultToFile(newProducts, gennum, fout)
        reactor.addCompounds(newProducts, reactions)
        results[gennum] = [Chem.MolToSmiles(mol) for mol in newProducts]
    if fout:
        fout.close()
