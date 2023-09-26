import math
import random
import os

ranseed=1234

class FakeModel:
    def __init__(self, NPars_Set, maxrank):
        self.NPars = NPars_Set
        self.randy = random.Random()
        self.randy.seed(ranseed)
        self.coefficient_sin = 50.0 * self.randy.random()
        self.coefficient_cos = 50.0 * self.randy.random()
        self.coefficient_exp = 50.0 * self.randy.random()
        self.LAMBDA = 2
        self.A = []
        self.Yname = ""
        self.Y = 0.0
        self.SigmaY = 1.0


    def GetY_1(self, iY, Yname, X):
        NPars = len(X)
        Lambda = 2.5
        arg = sum([self.randy.random() * x / (2.0 * math.pi) for x in X])
        Y = self.coefficient_sin * math.sin(arg / Lambda) + self.coefficient_cos * math.cos(arg / Lambda)
        arg = sum([self.randy.random() * x / (2.0 * math.pi) for x in X])
        Y += self.coefficient_exp * math.exp(arg / (Lambda * NPars))
        SigmaY = 1.0

        self.Yname = Yname
        self.Y = Y
        self.SigmaY = SigmaY
        return Y, SigmaY


with open("Info/modelpar_info.txt", "r") as file:
    NPars = sum(1 for line in file)

X = [0] * NPars
maxrank = 1 - NPars
SigmaReal = 10

Ynames = []
with open("Info/observable_info.txt", "r") as file1:
    for line_ob in file1:
        parts = line_ob.strip().split()
        if parts:  # Check if the line is not empty
            Yname = parts[0]
            Ynames.append(Yname)

fakeModel = FakeModel(NPars, maxrank)


itrain = 0
while True:
    filename = f"modelruns/run{itrain}/mod_parameters.txt"

    if not os.path.exists(filename):
        break
    print()
    with open(filename, 'r') as f:
        for ipar in range(NPars):
            line = f.readline().strip()
            if not line:
                continue
            try:
                parname, value = line.split()
                X[ipar] = float(value)
            except ValueError:
                print(f"Problematic line in file: {line}")
                continue

    filename = f"modelruns/run{itrain}/obs.txt"
    with open(filename, 'w') as f:
        for iY in range(len(Ynames)):
            Yname = Ynames[iY]
            Y, SigmaY = fakeModel.GetY_1(iY, Yname, X)
            print(f"Writing: {Yname} {Y} {SigmaY}")
            f.write(f"{Yname} {Y} {SigmaY}\n")
    itrain += 1
