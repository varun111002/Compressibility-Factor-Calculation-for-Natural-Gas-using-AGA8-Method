import aga8
import pprint
def example(t,p):
    #T = 273.15 + ((t - 32.0) * (5.0/9.0))
    if __name__ == '__main__':
        gasCompose = {
            'methane': 96.5222,
            'nitrogen': 0.2595,
            'carbonDioxide': 0.5956,
            'ethane': 1.8186,
            'propane': 0.4596,
            'water': 0,
            'hydrogenSulfide': 0,
            'hydrogen': 0,
            'carbonMonoxide': 0,
            'oxygen': 0,
            'iButane': 0.0977,
            'nButane': 0.1007,
            'iPentane': 0.0473,
            'nPentane': 0.0324,
            'nHexane': 0.0664,
            'nHeptane': 0.0,
            'nOctane': 0,
            'nNonane': 0,
            'nDecane': 0,
            'helium': 0.01,
            'argon': 0
            }
    a = aga8.AGA8()
    result = a.CalculateZ(gasCompositions=gasCompose, temperatureInKelvin=t, pressurePSI=p)
    return result
#set up temperature array
T= [32.00, 32.00, 32.00, 32.00, 32.00, 32.00, 32.00, 32.00]

#set up pressure array
P = [14.73, 100.00, 200.00, 400.00, 600.00, 800.00, 1000.00, 1200.00]
result = []
for i in T:
    if i is None:
        i = 0.0
        print(float(i))
    a = ((i - 32.0) * (5.0/9.0))+273.15
    for j in P:
        if j is None:
            j = 0.0
            print(float(j))
    print(i, 'F',a,'K',j, 'psia')
    pprint.pprint(example(a,j))
