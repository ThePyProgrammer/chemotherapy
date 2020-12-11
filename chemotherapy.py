# chemotherapy
prefixes = """meth
eth
prop
but
pent
hex
hept
oct
non
dec
undec
dodec
tridec
tetradec
pentadec
hexadec
heptadex
octadex
nonadec
icos
henicos
docos
tricos
tetracos
pentacos
hexacos
heptacos
octacos
nonocos
triacont
hentriacont
dotriacont
tritiacont
tetriacont
pentatriacont
hexatriacont
heptatriacont
octatriacont
nonatriacont
tetracont
hentetracont
dotetracont
tritetracont
tetratetracont
pentatetracont
hexatetracont
heptatetracont
octatetracont
nonatetracont
pentacont
henpentacont
dopentacont
tripentacont
tetrapentacont
pentapentacont
hexapentacont
heptapentacont
octapentacont
nonapentacont
hexacont
henhexacont
dohexacont
trihexacont
tetrahexacont
pentahexacont
hexahexacont
heptahexacont
octahexacont
nonahexacont
heptacont
henheptacont
doheptacont
triheptacont
tetraheptacont
pentaheptacont
hexaheptacont
heptaheptacont
octaheptacont
nonaheptacont
octacont
henoctacont
dooctacont
trioctacont
tetraoctacont
pentaoctacont
hexaoctacont
heptaoctacont
octaoctacont
nonaoctacont
nonacont
hennonacont
dononacont
trinonacont
tetranonacont
pentanonacont
hexanonacont
heptanonacont
octanonacont
nonanonacont
hect
henihect
dohect
trihect
tetrahect
pentahect
hexahect
heptahect
octahect
nonahect
decahect
undecahect
dodecahect
tridecahect
tetradecahect
pentadecahect
hexadecahect
heptadecahect
octadecahect
nonadecahect
icosahect""".split('\n')

class Atom:
    bd = 4
    def __init__(self):
        self.bondedto = [None for i in range(self.bd+1)]

    def attachBD(self, index, nextatom, bondtype=1):
        if index < self.bd:
            for i in range(bondtype):
                self.bondedto[index-1+bondtype] = nextatom
                if type(nextatom) == Hydrogen:
                    nextatom.attachBD(self)
        return self

    def attachLP(self, index, lp):
        if index < self.bd:
            self.bondedto[index] = lp
        return self

class Carbon(Atom):
    #bd = 4
    def __init__(self):
        super().__init__()
        for i in range(1, 3):
            self.bondedto[i] = LonePair()

class Hydrogen(Atom):
    bd = 1
    def __init__(self):
        super().__init__()
        
    def attachBD(self, nextatom):
        self.bondedto[1] = nextatom
        if type(nextatom) == Hydrogen:
            nextatom.attachBD(self)
        return self
        

class LonePair: pass

class Alkane:
    def __init__(self, n=None):
        if self.n == None:
            self.n = self.C = self.H = 'unknown'
        elif self.n == 1:
            self.n = self.C = n
            self.H = 2*n + 2
            self.mainC = Carbon()
            self.Hs = [Hydrogen() for i in range(self.H)]
            self.mainC.attachBD(1, self.Hs[0]).attachBD(2, self.Hs[1]).attachBD(3, self.Hs[2]).attachBD(4, self.Hs[3])
        else:
            self.n = self.C = n
            self.H = 2*n + 2
            self.Cs = [Carbon() for i in range(self.C)]
            self.Hs = [Hydrogen() for i in range(self.H)]
            self.Cs[0].attachBD(1, self.Hs[0]).attachBD(2, self.Hs[1]).attachBD(3, self.Hs[2]).attachBD(4, self.Cs[1])
            self.Cs[-1].attachBD(4, self.Hs[-1]).attachBD(3, self.Hs[-2]).attachBD(2, self.Hs[-3]).attachBD(1, self.Cs[-2])
            for i in range(1, len(self.Cs)-1):
                self.Cs[i].attachBD(3, self.Hs[i*2-1]).attachBD(2, self.Hs[i*2]).attachBD(1, self.Cs[i-1]).attachBD(4, self.Cs[i+1])
    def placeAt(self, Cindex, BDindex, atom, bondtype=1):
        self.Cs[Cindex].attachBD(BDindex, atom, bondtype)

    

class Alkene:
    def __init__(self, n=None, doubleind=1):
        if self.n == None:
            self.n = self.C = self.H = 'unknown'
        else:
            self.n = self.C = n
            self.H = 2*n
            self.Cs = [Carbon() for i in range(self.C)]
            self.Hs = [Hydrogen() for i in range(self.H)]
            if self.doubleind == 1:
                self.Cs[0].attachBD(1, self.Hs[0]).attachBD(2, self.Hs[1]).attachBD(3, self.Cs[1], 2)
                if n == 1:
                    self.Cs[-1].attachBD(4, self.Hs[-1]).attachBD(3, self.Hs[-2]).attachBD(1, self.Cs[-2], 2)
                elif n > 2:
                    self.Cs[1].attachBD(1, self.Cs[0], 2).attachBD(3, self.Hs[2]).attachBD(4, self.Cs[2])
                    self.Cs[-1].attachBD(4, self.Hs[-1]).attachBD(3, self.Hs[-2]).attachBD(2, self.Hs[-3]).attachBD(1, self.Cs[-2])
                    for i in range(2, len(self.Cs)-1):
                        self.Cs[i].attachBD(3, self.Hs[i*2-1]).attachBD(2, self.Hs[i*2]).attachBD(1, self.Cs[i-1]).attachBD(4, self.Cs[i+1])

            elif self.doubleind == n-1:
                self.Cs[-1].attachBD(4, self.Hs[-1]).attachBD(3, self.Hs[-2]).attachBD(1, self.Cs[-2], 2)
                if n == 1:
                    self.Cs[0].attachBD(1, self.Hs[0]).attachBD(2, self.Hs[1]).attachBD(3, self.Cs[1], 2)
                elif n > 2:
                    self.Cs[-2].attachBD(3, self.Cs[0], 2).attachBD(1, self.Cs[2]).attachBD(2, self.Hs[-3])
                    self.Cs[0].attachBD(1, self.Hs[0]).attachBD(2, self.Hs[1]).attachBD(3, self.Hs[2]).attachBD(4, self.Cs[1])
                    for i in range(2, len(self.Cs)-1):
                        self.Cs[i].attachBD(3, self.Hs[i*2-1]).attachBD(2, self.Hs[i*2]).attachBD(1, self.Cs[i-1]).attachBD(4, self.Cs[i+1])

    def placeAt(self, Cindex, BDindex, atom, bondtype=1):
        self.Cs[Cindex].attachBD(BDindex, atom, bondtype)
        


class Molecule:
    def __init__(self, formula):
        self.formula = formula.split('-')
        if formula[-1].strip().endswith("ane"):
            for i in range(len(prefixes)):
                if formula[-1].strip()[:-3] == prefixes[i].strip():
                    self.base = Alkane(formula)
            self.base = Alkane()
            
