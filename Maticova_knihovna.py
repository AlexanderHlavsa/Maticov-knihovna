import math #pouze pro zaokrouhlovani a odmocniny
import cmath
def zaokrouhleni(cislo):
    if math.isnan(cislo.real) or math.isnan(cislo.imag):
        return cislo  
    if math.isclose(abs(cislo.real), round(abs(cislo.real)), abs_tol=1e-6) and math.isclose(abs(cislo.imag), round(abs(cislo.imag)), abs_tol=1e-6):
        cislo = (round(cislo.real) + round(cislo.imag) * 1j) if cislo.imag!= 0 else round(cislo.real)
    return cislo

def round_complex(cislo, pocet):
    return complex(round(cislo.real, pocet),round(cislo.imag, pocet)) if round(cislo.imag, pocet) != 0 else round(cislo.real, pocet)

class Matrix:
    ###zaklad
    def __init__(self, data):
        self.data = data
        self.radky = len(data)
        self.sloupce = len(data[0]) if self.radky > 0 else 0
        self.prohozeni = 1 #determinant

    def __str__(self):
        copy = self.copy()
        for radek in range(copy.radky):
            for cislo in range(copy.sloupce):
                copy.data[radek][cislo] = round_complex(zaokrouhleni(copy.data[radek][cislo]),3)
        return "\n".join("  ".join(map(str, radek)) for radek in copy.data)
        
    def __repr__(self):
        copy = self.copy()
        for radek in range(copy.radky):
            for cislo in range(copy.sloupce):
                copy.data[radek][cislo] = round_complex(zaokrouhleni(copy.data[radek][cislo]),3)
        return f"Matrix({copy.data})"

    def copy(self):
        return Matrix([[self.data[i][j] for j in range(self.sloupce)] for i in range(self.radky)])
    
    ###


    #nacteni matice
    @staticmethod
    def from_list(list, m, n):
        if m*n != len(list):
            raise ValueError("Délka seznamu neodpovídá rozměrům matice.")
        k = 0
        matice = []
        for _ in range(m):
            mat = []
            for _ in range(n):
                mat.append(list[k])
                k +=1
            matice.append(mat)
        return Matrix(matice)
    
    def nul(self,radky,sloupce):
        data = [[0]*sloupce for _ in range(radky)]
        return Matrix(data)
    #typ matice
    def type(self):
        types = []
        if self.nulova(): 
            types.append("nulova")
            if self.ctvercova: types.append("ctvercova")
        elif self.jednotkova(): types.append("jednotkova, diagonalni, symetricka, ctvercova")
        elif self.diagonalni(): types.append("diagonalni, symetricka, ctvercova")
        elif self.symetricka(): types.append("symetricka, ctvercova")
        elif self.ctvercova(): types.append("ctvercova")
        if self.regularni(): types.append("regularni")
        return ', '.join(types)
        
    def ctvercova(self):
        if self.radky == self.sloupce:
            return True
        else: return False

    def diagonalni(self):
        if not self.ctvercova():
            return False
        for i in range (self.radky):
            for j in range(self.sloupce):
                if i != j and self.data[i][j] != 0:
                    return False
        return True
    def jednotkova(self):
        if not self.diagonalni():
            return False
        for i in range(self.radky):
            if self.data[i][i] != 1:
                return False
        return True
        
    def symetricka(self):
        if not self.ctvercova():
            return False
        for i in range(self.radky):
            for j in range(i, self.sloupce):
                if i != j and self.data[i][j] != self.data[j][i]:
                    return False
        return True
    
    def nulova(self):
        return all(all(abs(cislo) < 1e-10 for cislo in row) for row in self.data)
    
    def regularni(self):
        if not self.ctvercova():
            return False
        copy = self.copy()
        if copy.pocet_pivotu() != self.sloupce:
            return False
        return True

    def horni_trojuhelnikova(self):
        if not self.ctvercova():
            return False
        for j in range(self.sloupce):
            for i in range(j+1, self.radky):
                if self.data[i][j] != 0:
                    return False
        return True

    #transpozice
    def transpozice(self):
        nova_data = [[0]*self.radky for _ in range(self.sloupce)]
        for i in range(self.radky):
            for j in range(self.sloupce):
                nova_data [j][i] = self.data[i][j]
        return Matrix(nova_data)

    #scitani matic
    def stejny_typ_a_rozmery(self, other):
        if not isinstance(other, Matrix):
            raise TypeError("Objekty nejsou třídy 'Matrix'.")
        if self.radky != other.radky or self.sloupce != other.sloupce:
            raise ValueError("Matice nemají stejný typ")
        
    def __add__(self,other):
        self.stejny_typ_a_rozmery(other)
        nova_data = [[self.data[i][j] + other.data[i][j] for j in range(self.sloupce)] for i in range(self.radky)]
        return Matrix(nova_data)
    
    def __sub__(self, other):
        self.stejny_typ_a_rozmery(other)
        nova_data = [[self.data[i][j] - other.data[i][j] for j in range(self.sloupce)] for i in range(self.radky)]
        return Matrix(nova_data)
    
    def __neg__(self):
        nova_data = [[-self.data[i][j] for j in range(self.sloupce)] for i in range(self.radky)]
        return Matrix(nova_data)
    
    #nasobeni matic vcetne konstant
    def moznost_nasobeni(self,other):
        if not isinstance(other, Matrix):
            raise TypeError("Objekty nejsou třídy 'Matrix' nebo 'int' nebo 'float'.")
        if self.sloupce != other.radky:
            raise ValueError("Matice nelze násobit")
        
    def __mul__(self,other):
        if isinstance(other, (int,float,complex)):
            nova_data = [[0]*self.sloupce for _ in range(self.radky)]
            for i in range(self.radky):
                for j in range(self.sloupce):
                    nova_data[i][j] = other*self.data[i][j]
            return Matrix(nova_data)
        else:
            self.moznost_nasobeni(other)
            nova_data = [[0]*other.sloupce for _ in range(self.radky)]
            for i in range(self.radky):
                for j in range(other.sloupce):
                    for k in range(self.sloupce):
                        nova_data[i][j] += self.data[i][k] * other.data[k][j]
            return Matrix(nova_data)
        
    def __rmul__(self,other):
        return self.__mul__(other)
    
    #kroneckeruv soucin
    def kron(self, other):
        nova_data = []
        for ai in range(self.radky):
            for bi in range(other.radky):
                radek = []
                for aj in range(self.sloupce):
                   for bj in range(other.sloupce):
                       radek.append(self.data[ai][aj]* other.data[bi][bj])
                nova_data.append(radek)
        return Matrix(nova_data)
    
    #Gaussova eliminace
    def Gauss(self):
        radek = 0
        prohozeni = 1
        for sloupec in range(self.sloupce):
            radek,prohozeni = self.eliminace_sloupce(sloupec, radek, prohozeni)
        self.prohozeni = prohozeni
        return self

    def eliminace_sloupce(self, sloupec, radek, prohozeni):
        data = self.data
        puv_radek = radek
        while radek < self.radky and data[radek][sloupec] == 0:
            radek += 1
        if radek == self.radky:
            return puv_radek, prohozeni
        if puv_radek != radek:
            data[puv_radek], data[radek] = data[radek], data[puv_radek]
            prohozeni *= -1 

        for i in range(puv_radek +1, self.radky):          
            if data[i][sloupec] != 0:
                nasobek = data[puv_radek][sloupec]/data[i][sloupec]
                if nasobek == 0: continue
                prohozeni *=1/nasobek
                for j in range(sloupec, self.sloupce):
                    data[i][j] = zaokrouhleni((data[i][j]*nasobek) - data[puv_radek][j])
        return puv_radek +1, prohozeni

    def Gauss_Jordan(self):
        self.Gauss()
        self.Jordan()
        return self
    
    def Jordan(self):
        data = self.data
        for radek in range(self.radky):
            sloupec = radek
            while sloupec < self.sloupce and data[radek][sloupec] == 0:
                sloupec += 1
            if sloupec >= self.sloupce:
                continue
            delitel = data[radek][sloupec]
            for j in range(sloupec, self.sloupce):
                data[radek][j] = data[radek][j]/delitel

            for i in range(radek):
                nasobek = data[i][sloupec]
                for j in range(sloupec, self.sloupce):
                    data[i][j] -= nasobek * data[radek][j]

            for i in range(radek + 1, self.radky):
                nasobek = data[i][sloupec]
                for j in range(sloupec, self.sloupce):
                    data[i][j] -= nasobek * data[radek][j]

        return self

    #pocitani s matici
    def pocet_pivotu(self): #= hodnost = pocet nenulovych radku = pocet bazovych sloupcu
        self.Gauss()
        data = self.data
        pocet_nul_radku = 0
        for radek in range(self.radky -1, -1, -1):
            for sloupec in range (self.sloupce):
                if data[radek][sloupec] != 0:
                    return self.radky - pocet_nul_radku
            pocet_nul_radku += 1
        return 0
    
    def pozice_pivotu(self):
        pocet_pivotu = self.pocet_pivotu()
        pozice = []
        sloupec = 0
        data = self.data    
        for i in range(self.radky):
            for j in range(sloupec, self.sloupce):
                if data[i][j] != 0:
                    pozice.append([i,j])
                    sloupec = j+1
                    break
            if len(pozice) == pocet_pivotu:
                break
        return pozice
    
    #hodnost
    def hodnost(self):
        copy = self.copy()
        return copy.pocet_pivotu()
    
    def dimenze_jadra(self):
        return self.sloupce - self.hodnost()

    def pridat_pravou_stranu(self, prava_strana):
        data = self.data
        for i in range(self.radky):
            for j in range(len(prava_strana.data[0])):
                data[i].append(prava_strana.data[i][j])
        self.radky = len(data)
        self.sloupce = len(data[0]) if self.radky > 0 else 0
        return 
    
    def reseni(self,prava_strana):
        copy = self.copy()
        if not isinstance(prava_strana, Matrix):
            raise TypeError("Pravá strana není třídy 'Matrix'.")
        if prava_strana.sloupce != 1 or prava_strana.radky != copy.radky:
            raise ValueError("Pravá strana nemá platný tvar.")
        radky = copy.radky
        sloupce = copy.sloupce
        copy.pridat_pravou_stranu(prava_strana)
        pozice_pivotu = copy.pozice_pivotu()
        copy.Jordan()
        bazove_sloupce = set()
        
        for i in range(len(pozice_pivotu)):
            bazove_sloupce.add(pozice_pivotu[i][1])
        if sloupce in bazove_sloupce:
            raise ValueError("Soustava nemá řešení")
        neb_sloupce = set(range(0,copy.sloupce))- bazove_sloupce
        neb_sloupce -= {copy.sloupce-1}

        vektor = [[False] for _ in range(self.sloupce)]
        for i in neb_sloupce:
            vektor[i][0] = 1
        for pivot in reversed(pozice_pivotu):
            radek = pivot[0]
            sloupec = pivot[1]
            vysledek = copy.data[radek][-1] - sum(copy.data[radek][j]*vektor[j][0] for j in range(sloupec+1,copy.sloupce-1))
            vektor[sloupec][0] = vysledek
        return Matrix(vektor)

    #jadro matice
    def Ker(self):
        copy = self.copy()
        prava_strana = Matrix([[0]for _ in range(self.radky)])
        copy.pridat_pravou_stranu(prava_strana)
        pozice_pivotu = copy.pozice_pivotu()
        copy.Jordan()
        bazove_sloupce = set()

        for i in range(len(pozice_pivotu)):
            bazove_sloupce.add(pozice_pivotu[i][1])
        sloupce = set(range(0,copy.sloupce))
        neb_sloupce = sloupce - bazove_sloupce
        neb_sloupce -= {copy.sloupce-1}
        Jadro = []
        for i in neb_sloupce:
            vektor = [0]*self.sloupce
            vektor[i] = 1
            for i in bazove_sloupce:
                vektor[i] = False
            for pivot in reversed(pozice_pivotu):
                radek = pivot[0]
                sloupec = pivot[1]
                vysledek = copy.data[radek][-1] - sum(copy.data[radek][j]*vektor[j] for j in range(sloupec+1,copy.sloupce-1))
                vektor[sloupec] = vysledek
            Jadro.append(vektor)
        for i in range(len(Jadro)):
            Jadro[i] = Matrix([[Jadro[i][j]]for j in range(self.sloupce)])
        if len(Jadro) == 0:
            return [Matrix([[0]for _ in range(self.sloupce)])]
        return Jadro

    def determinant(self):
        copy = self.copy()
        if not copy.ctvercova():
            raise ValueError ("Matice není čtvercová, nelze určit determinant")
        copy.Gauss()
        determinant = copy.prohozeni
        for i in range(copy.radky):
            determinant *= copy.data[i][i]
        return round_complex(zaokrouhleni(determinant),3)
    
    def jednotk(self):
        data = [[0]*self.radky for _ in range(self.radky)]
        for i in range(self.radky):
            data[i][i] = 1
        return Matrix(data)

    def inverz(self):
        if not self.regularni():
            raise ValueError ("Matice není regulární a nemá inverz")
        copy = self.copy()
        copy.pridat_pravou_stranu(copy.jednotk())
        copy = Matrix(copy.data)
        copy.Gauss_Jordan()
        nova_data = [[0]*copy.radky for _ in range(copy.radky)]
        for i in range(copy.radky):
            for j in range(copy.radky, copy.sloupce):
                nova_data[i][j-copy.radky] = copy.data[i][j]

        return Matrix(nova_data)
    ### QR rozklad
    def skalarni_soucin(self,v2):
        if not isinstance(v2, Matrix):
            raise TypeError("Objekty nejsou třídy 'Matrix'.")
        if self.sloupce != 1 and v2.sloupce != 1:
            raise ValueError("První i druhý objekt musí být vektory (musí mít jeden sloupec)")
        elif self.sloupce != 1:
            raise ValueError("První objekt musí být vektor (musí mít jeden sloupec)")
        elif v2.sloupce != 1:
            raise ValueError("První objekt musí být vektor (musí mít jeden sloupec)")
        elif v2.radky != self.radky:
            raise ValueError("Vektory musí mít stejný počet řádků")
        return sum(self.data[i][0].conjugate() * v2.data[i][0] for i in range(self.radky))

    
    def LN_sloupce(self):
        copy = self.copy()
        if copy.sloupce != copy.pocet_pivotu():
            return False
        return True
    
    def norma(self):
        return math.sqrt(sum(abs(self.data[i][0])**2 for i in range(self.radky)))


    def normalizace(self):
        norma = self.norma()
        if norma != 0 and norma != 1:
            for i in range(self.radky):
                self.data[i][0] /= norma
        return self
    
    def ortogonalni(self):
        for j in range(self.sloupce):
            sloupec = Matrix([[self.data[i][j]] for i in range(self.radky)])
            if zaokrouhleni(sloupec.skalarni_soucin(sloupec)) != 1:
                return False
        return True

    ### A = Q*R
    def QR(self,stav = True):
        if stav:
            if not self.LN_sloupce():
                raise ValueError("Matice nemá lineárně nezávislé sloupce")
        copy = self.copy()
        Qdata = [[0]*self.sloupce for _ in range(self.radky)]
        Rdata = [[0]*self.sloupce for _ in range(self.sloupce)]
        prvni_vektor = Matrix([[copy.data[i][0]] for i in range(copy.radky)])
        norma = prvni_vektor.norma()
        prvni_vektor *= (1/norma) if norma != 0 else 0
        for i in range(self.radky):
            Qdata[i][0] = prvni_vektor.data[i][0]
        Rdata[0][0] = norma
        for vektor in range(1, self.sloupce):
            w = Matrix([[0] for _ in range(copy.radky)])
            v = Matrix([[copy.data[s][vektor]] for s in range(self.radky)])
            for j in range(vektor):
                uj = Matrix([[Qdata[s][j]] for s in range(self.radky)])
                sk_soucin = v.skalarni_soucin(uj)
                w += (uj * sk_soucin)
                Rdata[j][vektor] = sk_soucin
            ui = v - w
            norma = ui.norma()
            ui *= (1/norma) if norma != 0 else 0
            Rdata[vektor][vektor] = norma
            for i in range(self.radky):
                Qdata[i][vektor] = ui.data[i][0]    
        Q = Matrix(Qdata)
        R = Matrix(Rdata)
        return Q, R
    #Cholevskeho rozklad
    def Forbien_norma(self):
        return math.sqrt(sum(sum(abs(prvek)**2 for prvek in row) for row in self.data))

    def vlastni_cisla(self, tol=1e-10, iterace=1000):
        if not self.ctvercova():
            raise ValueError("Matice není čtvercová")
        copy = self.copy()
        vlastni_cisla = []
        for _ in range(iterace):
            Q, R = copy.QR(False)
            copy = R * Q
            
            off_diagonal_norm = (copy-(Matrix([[copy.data[i][i] if i == j else 0 for j in range(self.sloupce)] for i in range(self.radky)]))).Forbien_norma()
            if off_diagonal_norm < tol:
                break
        for cislo in [zaokrouhleni(copy.data[i][i]) for i in range(self.radky)]:
            if self.overeni(cislo) is not None:
                vlastni_cisla.append(cislo)
        return vlastni_cisla
    
    def overeni(self,cislo):
        for vektor in (self - cislo*self.jednotk()).Ker():
            if vektor.nulova() or (self*vektor).zaokr_matice().data != (cislo*vektor).zaokr_matice().data:
                return None
        return cislo

    def vlastni_vektory(self):
        vlastni_vektory = []
        for v_cislo in set(self.vlastni_cisla()):
            for vektor in (self - v_cislo*self.jednotk()).Ker():
                vlastni_vektory.append(vektor)
        return vlastni_vektory
    
    def zaokr_matice(self):
        copy = self.copy()
        for radek in range(copy.radky):
            for cislo in range(copy.sloupce):
                copy.data[radek][cislo] = round_complex(zaokrouhleni(copy.data[radek][cislo]),4)
        return copy

    def pozitivne_definitni(self):
        if not self.symetricka():
            return False
        for vlastni_cislo in self.vlastni_cisla():
            if vlastni_cislo <= 0:
                return False
        return True
    
    def Choleskeho_rozklad(self):
        if not self.pozitivne_definitni():
            raise ValueError("Matice není pozitivně definitní")
        Ldata = [[0]*self.sloupce for _ in range(self.radky)]
        for radek in range(self.radky):
            for sloupec in range(radek + 1):
                sum_k = sum(Ldata[radek][k] * Ldata[sloupec][k] for k in range(sloupec))
                if sloupec == radek:
                    Ldata[radek][sloupec] = (self.data[radek][radek] - sum_k) ** 0.5
                else:
                    Ldata[radek][sloupec] = (self.data[radek][sloupec] - sum_k) / Ldata[sloupec][sloupec]
        return Matrix(Ldata)
