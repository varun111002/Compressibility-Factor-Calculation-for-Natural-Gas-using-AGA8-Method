import pprint

import aga_const
import math

gas_list = ['methane', 'nitrogen', 'carbonDioxide', 'ethane', 'propane', 'water', 'hydrogenSulfide', 'hydrogen',
            'carbonMonoxide', 'oxygen', 'iButane', 'nButane', 'iPentane', 'nPentane', 'nHexane', 'nHeptane', 'nOctane',
            'nNonane', 'nDecane', 'helium', 'argon']


class AGA8:
    RGAS = 0.00831451

    EPSR = EPSMIN = RHO = RHOL = RHOH = PRHOL = PRHOH = 0

    FN = []

    UU = None
    RK3PO = None
    Q2PO = None
    BMIX = None

    CODE = None

    RK5PO = None
    RK2P5 = None
    U5PO = None
    U2P5 = None
    WW = None
    Q1PO = None
    HH = None

    MWX = None

    B = []

    TOLD = None

    gasComps = []

    def __init__(self):
        pass

    def CalculateZ(self, gasCompositions, temperatureInKelvin, pressurePSI):

        if gasCompositions is None:
            raise ValueError("gasComps must have a value")

        # If dict value, convert it to a list
        self.gasComps = []
        if type(gasCompositions) == dict:
            try:
                for ii in range(len(gas_list)):
                    self.gasComps.append(gasCompositions[gas_list[ii]] / 100.0)
            except KeyError as e:
                print (e)
                return None
        elif type(gasCompositions) != list:
            raise ValueError('gasComps must be a list or dictionary')
        else:
            self.gasComps = [g / 100.0 for g in gasCompositions]

        # Check total percentage.
        if not (.99 < sum(self.gasComps) < 1.01):
            raise ValueError("Total Gas Composition Does not add up to 100%: {}".format(sum(self.gasComps)))

        self.B = [0.0] * 18
        self.RK5PO = 0.0
        self.RK2P5 = 0.0
        self.U5PO = 0.0
        self.U2P5 = 0.0
        self.WW = 0.0
        self.Q1PO = 0.0
        self.HH = 0.0
        self.MWX = 0
        for i in range(len(self.gasComps)):
            self.MWX += (self.gasComps[i] * aga_const.MolarMass[i])

        for i in range(len(self.gasComps)):
            self.RK2P5 += float(self.gasComps[i]) * aga_const.K_Size[i] * aga_const.K_Size[i] * \
                          math.sqrt(aga_const.K_Size[i])
            self.U2P5 += float(self.gasComps[i]) * aga_const.E_Energy[i] * aga_const.E_Energy[i] * math.sqrt(
                aga_const.E_Energy[i])
            self.HH += float(self.gasComps[i]) * aga_const.W_Association[i]
            self.Q1PO += float(self.gasComps[i]) * aga_const.Q_Quadrupole[i]
            self.WW += float(self.gasComps[i]) * aga_const.G_Orientation[i]
            for j in range(i, len(self.gasComps)):
                if i != j:
                    XIJ = 2.0 * self.gasComps[i] * self.gasComps[j]
                else:
                    XIJ = self.gasComps[i] * self.gasComps[j]

                if aga_const.K_BinaryInteration[i][j] != 1:
                    self.RK5PO += XIJ * (math.pow(aga_const.K_BinaryInteration[i][j], 5) - 1) * math.sqrt(
                        math.pow(aga_const.K_Size[i], 5) * math.pow(aga_const.K_Size[j], 5))

                if aga_const.U_BinaryInteration[i][j] != 1:
                    self.U5PO += XIJ * (math.pow(aga_const.U_BinaryInteration[i][j], 5) - 1) * math.sqrt(
                        (math.pow(aga_const.E_Energy[i], 5) * math.pow(aga_const.E_Energy[j], 5)))

                if aga_const.G_BinaryInteration[i][j] != 1:
                    self.WW += XIJ * (aga_const.G_BinaryInteration[i][j] - 1) * \
                               ((aga_const.G_Orientation[i] + aga_const.G_Orientation[j]) / 2.0)

                EIJ = aga_const.E_BinaryInteration[i][j] * math.sqrt(aga_const.E_Energy[i] * aga_const.E_Energy[j])
                WIJ = aga_const.G_BinaryInteration[i][j] * \
                      (aga_const.G_Orientation[i] + aga_const.G_Orientation[j]) / 2.0
                EOP5 = math.sqrt(EIJ)
                E2PO = EIJ * EIJ
                E3PO = EIJ * E2PO
                E3P5 = E3PO * EOP5
                E4P5 = EIJ * E3P5
                E6PO = E3PO * E3PO
                E11PO = E4P5 * E4P5 * E2PO
                E7P5 = E4P5 * EIJ * E2PO
                E9P5 = E7P5 * E2PO
                E12PO = E11PO * EIJ
                E12P5 = E12PO * EOP5
                S3 = XIJ * math.sqrt(math.pow(aga_const.K_Size[i], 3) * math.pow(aga_const.K_Size[j], 3))
                self.B[0] += S3
                self.B[1] += S3 * EOP5
                self.B[2] += S3 * EIJ
                self.B[3] += S3 * E3P5
                self.B[4] += S3 * WIJ / EOP5
                self.B[5] += S3 * WIJ * E4P5
                self.B[6] += S3 * aga_const.Q_Quadrupole[i] * aga_const.Q_Quadrupole[j] * EOP5
                self.B[7] += S3 * aga_const.S_Dipole[i] * aga_const.S_Dipole[j] * E7P5
                self.B[8] += S3 * aga_const.S_Dipole[i] * aga_const.S_Dipole[j] * E9P5
                self.B[9] += S3 * aga_const.W_Association[i] * aga_const.W_Association[j] * E6PO
                self.B[10] += S3 * aga_const.W_Association[i] * aga_const.W_Association[j] * E12PO
                self.B[11] += S3 * aga_const.W_Association[i] * aga_const.W_Association[j] * E12P5
                self.B[12] += S3 * aga_const.F_HighTemp[i] * aga_const.F_HighTemp[j] / E6PO
                self.B[13] += S3 * E2PO
                self.B[14] += S3 * E3PO
                self.B[15] += S3 * aga_const.Q_Quadrupole[i] * aga_const.Q_Quadrupole[j] * E2PO
                self.B[16] += S3 * E2PO
                self.B[17] += S3 * E11PO

        for i in range(len(self.B)):
            self.B[i] = self.B[i] * aga_const.a_EquationOfState[i]

        self.RK3PO = math.pow((self.RK5PO + self.RK2P5 * self.RK2P5), 0.6)
        self.UU = math.pow((self.U5PO + self.U2P5 * self.U2P5), 0.2)
        self.Q2PO = self.Q1PO * self.Q1PO

        self.Temp(temperatureInKelvin)

        PF = pressurePSI * 6894.757 / 1000000
        DF = self.DDetail(PF, temperatureInKelvin)
        ZF = self.ZDetail(DF, temperatureInKelvin)

        aga8CalcResult = dict()

        aga8CalcResult['FlowCompressiblity'] = round(ZF, 10)

        return ZF

    def molarMass(self, gasComps):
        mass = 0
        for i in range(len(gasComps)):
            mass = mass + (aga_const.MolarMass[i] * gasComps[i])
        return mass

    def DDetail(self, Pressure, TempInKevin):

        IMAX = 150

        self.RHO = 0
        self.RHOL = 0
        self.RHOH = 0
        self.PRHOL = 0
        self.PRHOH = 0

        EPSP = 0.000001
        EPSR = 0.000001
        EPSMIN = 0.0000001

        self.BRACKET(TempInKevin, Pressure)
        if self.CODE in [1, 3]:
            return self.RHO

        X1 = self.RHOL
        X2 = self.RHOH
        Y1 = self.PRHOL - Pressure
        Y2 = self.PRHOH - Pressure
        DELX = X1 - X2
        DELPRV = DELX

        X3 = X1
        Y3 = Y1

        for i in range(IMAX):
            if (Y2 * Y3) > 0:
                X3 = X1
                Y3 = Y1
                DELX = X1 - X2
                DELPRV = DELX

            if math.fabs(Y3) < math.fabs(Y2):
                X1 = X2
                X2 = X3
                X3 = X1
                Y1 = Y2
                Y2 = Y3
                Y3 = Y1
            DELMIN = EPSMIN * math.fabs(X2)

            DELBIS = 0.5 * (X3 - X2)

            if math.fabs(DELPRV) < DELMIN or math.fabs(Y1) < math.fabs(Y2):
                DELX = DELBIS
                DELPRV = DELBIS
            else:
                if X3 != X1:
                    Y2MY3 = Y2 - Y3
                    Y3MY1 = Y3 - Y1
                    Y1MY2 = Y1 - Y2
                    XDENOM = -Y1MY2 * Y2MY3 * Y3MY1
                    XNUMBER = X1 * Y2 * Y3 * Y2MY3 + X2 * Y3 * Y1 * Y3MY1 + X3 * Y1 * Y2 * Y1MY2 - X2 * XDENOM
                else:
                    XNUMBER = (X2 - X1) * Y2
                    XDENOM = Y1 - Y2

                if 2 * math.fabs(XNUMBER) < math.fabs(DELPRV * XDENOM):
                    DELPRV = DELX
                    DELX = XNUMBER / XDENOM
                else:
                    DELX = DELBIS
                    DELPRV = DELBIS

                if (math.fabs(Y2) < EPSP * Pressure) and math.fabs(DELX) < EPSR * math.fabs(X2):
                    return X2 + DELX

            if math.fabs(DELX) < DELMIN:
                SGNDEL = DELBIS / math.fabs(DELBIS)
                DELX = 1.0000009 * SGNDEL * DELMIN
                DELPRV = DELX

            BOUNDN = DELX * (X2 + DELX - X3)

            if BOUNDN > 0:
                DELX = DELBIS
                DELPRV = DELBIS

            X1 = X2
            Y1 = Y2
            X2 += DELX
            Y2 = self.PDetail(X2, TempInKevin) - Pressure

        # If we exceed max number, then return X2
        return X2

    def BRACKET(self, TemperatureInKelvin, Pressure):
        self.CODE = 0
        IMAX = 200

        self.RHOL = 0
        self.PRHOL = 0
        self.RHOH = 0
        self.PRHOH = 0
        self.RHO = 0

        RHO1 = 0
        P1 = 0
        RHOMAX = 1.0 / self.RK3PO

        # if TemperatureInKelvin == 330.15:
        #     fdslkj = 2

        if TemperatureInKelvin > 1.2593 * self.UU:
            RHOMAX *= 20.0

        VIDEAL = self.RGAS * TemperatureInKelvin / Pressure

        if math.fabs(self.BMIX) < (0.167 * VIDEAL):
            RHO2 = 0.95 / (VIDEAL + self.BMIX)
        else:
            RHO2 = 1.15 / VIDEAL

        DEL = RHO2 / 20.0
        IT = 0

        while True:
            IT += 1
            if IT > IMAX:
                self.CODE = 3
                self.RHO = RHO2
                return

            if self.CODE != 2 and RHO2 > RHOMAX:
                self.CODE = 2
                DEL = 0.01 * (RHOMAX - RHO1) + Pressure / (self.RGAS * TemperatureInKelvin) / 20.0
                RHO2 = RHO1 + DEL
            else:
                P2 = self.PDetail(RHO2, TemperatureInKelvin)

                if P2 > Pressure:
                    self.RHOL = RHO1
                    self.PRHOL = P1
                    self.RHOH = RHO2
                    self.PRHOH = P2
                    return
                elif P2 > P1 and self.CODE == 2:
                    RHO1 = RHO2
                    P1 = P2
                    RHO2 = RHO1 + DEL
                elif P2 > P1 and self.CODE == 0:
                    DEL *= 2
                    RHO1 = RHO2
                    P1 = P2
                    RHO2 = RHO1 + DEL
                else:
                    self.CODE = 1
                    self.RHO = RHO1
                    return

    def PDetail(self, Rho, temp):
        return self.ZDetail(Rho, temp) * Rho * self.RGAS * temp

    def ZDetail(self, density, temperatureInKelvin):

        if temperatureInKelvin != self.TOLD:
            self.Temp(temperatureInKelvin)

        D1 = self.RK3PO * density
        D2 = D1 * D1
        D3 = D2 * D1
        D4 = D3 * D1
        D5 = D4 * D1
        D6 = D5 * D1
        D7 = D6 * D1
        D8 = D7 * D1
        D9 = D8 * D1

        EXP1 = math.exp(-D1)
        EXP2 = math.exp(-D2)
        EXP3 = math.exp(-D3)
        EXP4 = math.exp(-D4)

        result = 1 + self.BMIX * density + self.FN[12] * D1 * (EXP3 - 1 - 3 * D3 * EXP3) + \
                 (self.FN[13] + self.FN[14] + self.FN[15]) * D1 * (EXP2 - 1 - 2 * D2 * EXP2) + \
                 (self.FN[16] + self.FN[17]) * D1 * (EXP4 - 1 - 4 * D4 * EXP4) + \
                 (self.FN[18] + self.FN[19]) * D2 * 2 + \
                 (self.FN[20] + self.FN[21] + self.FN[22]) * D2 * (2 - 2 * D2) * EXP2 + \
                 (self.FN[23] + self.FN[24] + self.FN[25]) * D2 * (2 - 4 * D4) * EXP4 + \
                 self.FN[26] * D2 * (2 - 4 * D4) * EXP4 + \
                 self.FN[27] * D3 * 3 + \
                 (self.FN[28] + self.FN[29]) * D3 * (3 - D1) * EXP1 + \
                 (self.FN[30] + self.FN[31]) * D3 * (3 - 2 * D2) * EXP2 + \
                 (self.FN[32] + self.FN[33]) * D3 * (3 - 3 * D3) * EXP3 + \
                 (self.FN[34] + self.FN[35] + self.FN[36]) * D3 * (3 - 4 * D4) * EXP4 + \
                 (self.FN[37] + self.FN[38]) * D4 * 4 + \
                 (self.FN[39] + self.FN[40] + self.FN[41]) * D4 * (4 - 2 * D2) * EXP2 + \
                 (self.FN[42] + self.FN[43]) * D4 * (4 - 4 * D4) * EXP4 + \
                 self.FN[44] * D5 * 5 + \
                 (self.FN[45] + self.FN[46]) * D5 * (5 - 2 * D2) * EXP2 + \
                 (self.FN[47] + self.FN[48]) * D5 * (5 - 4 * D4) * EXP4 + \
                 self.FN[49] * D6 * 6 + \
                 self.FN[50] * D6 * (6 - 2 * D2) * EXP2 + \
                 self.FN[51] * D7 * 7 + \
                 self.FN[52] * D7 * (7 - 2 * D2) * EXP2 + \
                 self.FN[53] * D8 * (8 - D1) * EXP1 + \
                 (self.FN[54] + self.FN[55]) * D8 * (8 - 2 * D2) * EXP2 + \
                 (self.FN[56] + self.FN[57]) * D9 * (9 - 2 * D2) * EXP2

        return result

    def BMix(self, tempInKelvin):

        TOP5 = math.sqrt(tempInKelvin)
        T2PO = tempInKelvin * tempInKelvin
        T3PO = tempInKelvin * T2PO
        T3P5 = T3PO * TOP5
        T4P5 = tempInKelvin * T3P5
        T6PO = T3PO * T3PO
        T11PO = T4P5 * T4P5 * T2PO
        T7P5 = T6PO * tempInKelvin * TOP5
        T9P5 = T7P5 * T2PO
        T12PO = T9P5 * TOP5 * T2PO
        T12P5 = T12PO * TOP5

        self.BMIX = self.B[0] + self.B[1] / TOP5 + self.B[2] / tempInKelvin + self.B[3] / T3P5 + self.B[4] * TOP5 + \
                    self.B[5] / T4P5 + self.B[6] / TOP5 + self.B[7] / T7P5 + self.B[8] / T9P5 + self.B[9] / T6PO + \
                    self.B[10] / T12PO + self.B[11] / T12P5 + self.B[12] * T6PO + self.B[13] / T2PO + \
                    self.B[14] / T3PO + self.B[15] / T2PO + self.B[16] / T2PO + self.B[17] / T11PO
        return self.BMIX

    def Temp(self, temperatureInKelvin):
        self.FN = [0] * 58
        self.BMix(temperatureInKelvin)

        TR = temperatureInKelvin / self.UU
        TROP5 = math.sqrt(TR)
        TR1P5 = TR * TROP5
        TR2PO = TR * TR
        TR3PO = TR * TR2PO
        TR4PO = TR * TR3PO
        TR5PO = TR * TR4PO
        TR6PO = TR * TR5PO
        TR7PO = TR * TR6PO
        TR8PO = TR * TR7PO
        TR9PO = TR * TR8PO
        TR11PO = TR6PO * TR5PO
        TR13PO = TR6PO * TR7PO
        TR21PO = TR9PO * TR9PO * TR3PO
        TR22PO = TR * TR21PO
        TR23PO = TR * TR22PO

        self.FN[12] = aga_const.a_EquationOfState[12] * self.HH * TR6PO
        self.FN[13] = aga_const.a_EquationOfState[13] / TR2PO
        self.FN[14] = aga_const.a_EquationOfState[14] / TR3PO
        self.FN[15] = aga_const.a_EquationOfState[15] * self.Q2PO / TR2PO
        self.FN[16] = aga_const.a_EquationOfState[16] / TR2PO
        self.FN[17] = aga_const.a_EquationOfState[17] / TR11PO
        self.FN[18] = aga_const.a_EquationOfState[18] * TROP5
        self.FN[19] = aga_const.a_EquationOfState[19] / TROP5
        self.FN[20] = aga_const.a_EquationOfState[20]
        self.FN[21] = aga_const.a_EquationOfState[21] / TR4PO
        self.FN[22] = aga_const.a_EquationOfState[22] / TR6PO
        self.FN[23] = aga_const.a_EquationOfState[23] / TR21PO
        self.FN[24] = aga_const.a_EquationOfState[24] * self.WW / TR23PO
        self.FN[25] = aga_const.a_EquationOfState[25] * self.Q2PO / TR22PO
        self.FN[26] = aga_const.a_EquationOfState[26] * self.HH * TR
        self.FN[27] = aga_const.a_EquationOfState[27] * self.Q2PO * TROP5
        self.FN[28] = aga_const.a_EquationOfState[28] * self.WW / TR7PO
        self.FN[29] = aga_const.a_EquationOfState[29] * self.HH * TR
        self.FN[30] = aga_const.a_EquationOfState[30] / TR6PO
        self.FN[31] = aga_const.a_EquationOfState[31] * self.WW / TR4PO
        self.FN[32] = aga_const.a_EquationOfState[32] * self.WW / TR
        self.FN[33] = aga_const.a_EquationOfState[33] * self.WW / TR9PO
        self.FN[34] = aga_const.a_EquationOfState[34] * self.HH * TR13PO
        self.FN[35] = aga_const.a_EquationOfState[35] / TR21PO
        self.FN[36] = aga_const.a_EquationOfState[36] * self.Q2PO / TR21PO
        self.FN[37] = aga_const.a_EquationOfState[37] * TROP5
        self.FN[38] = aga_const.a_EquationOfState[38]
        self.FN[39] = aga_const.a_EquationOfState[39] / TR2PO
        self.FN[40] = aga_const.a_EquationOfState[40] / TR7PO
        self.FN[41] = aga_const.a_EquationOfState[41] * self.Q2PO / TR9PO
        self.FN[42] = aga_const.a_EquationOfState[42] / TR22PO
        self.FN[43] = aga_const.a_EquationOfState[43] / TR23PO
        self.FN[44] = aga_const.a_EquationOfState[44] / TR
        self.FN[45] = aga_const.a_EquationOfState[45] / TR9PO
        self.FN[46] = aga_const.a_EquationOfState[46] * self.Q2PO / TR3PO
        self.FN[47] = aga_const.a_EquationOfState[47] / TR8PO
        self.FN[48] = aga_const.a_EquationOfState[48] * self.Q2PO / TR23PO
        self.FN[49] = aga_const.a_EquationOfState[49] / TR1P5
        self.FN[50] = aga_const.a_EquationOfState[50] * self.WW / TR5PO
        self.FN[51] = aga_const.a_EquationOfState[51] * self.Q2PO * TROP5
        self.FN[52] = aga_const.a_EquationOfState[52] / TR4PO
        self.FN[53] = aga_const.a_EquationOfState[53] * self.WW / TR7PO
        self.FN[54] = aga_const.a_EquationOfState[54] / TR3PO
        self.FN[55] = aga_const.a_EquationOfState[55] * self.WW
        self.FN[56] = aga_const.a_EquationOfState[56] / TR
        self.FN[57] = aga_const.a_EquationOfState[57] * self.Q2PO

        return 0
