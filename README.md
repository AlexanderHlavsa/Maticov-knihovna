# Maticová knihovna
Alexander Hlavsa - Zápočtový projekt programování 2

Návod na instalaci
-	K instalaci maticové knihovny je třeba stáhnout soubor s kódem např. pod názvem Maticova_knihovna.py
-	Vložit tento soubor do stejné složky jako nový kód, který se chystáte psát
-	Otevřít tuto složku např. ve VS Code 
-	Do  nového kódu napsat příkaz „import Maticova_knihovna“ nebo „ import Maticova_knihovna as M“
-	Jakmile kód spustíte, VS code automaticky vytvoří složku __pycache__ a váš kód bude fungovat
-	Můžete začít psát

Popis knihovny
-	Základem celé knihovny je class Matrix, ve kterém jsou implementovány jednotlivé funkce
-	Všechny funkce lze volat buď A.funkce(), nebo mp.Matrix.funkce(A), kde A je matice pro kterou chceme funkci provést
-	Pro všechny funkce je zřízené vhodné zaokrouhlení pro pohodlnou a přehlednou práci: 
- Kód předchází pythonovským zaokrouhlovacím chybám tedy např. číslo 1.000000000001 se uloží a vypíše jako 1 
-	Při výpisu matice a jejích hodnot se vypíšou hodnoty zaokrouhlené na 3 desetinná místa, ale uloženy zůstávají přesné nezaokrouhlené hodnoty pro další přesné počítání (funkce __str__  a  __repr__)
- Každá matice A má 4 základní artributy:
  
  A.data – hodnoty matice v podobě 2D seznamu
  
  A.radky – počet řádků matice
  
  A.sloupce – počet sloupců matice
  
  A.prohozeni – slouží pouze pro počítání determinantu a její použití pro uživatele nemá roli
  
-	Knihovna obsahuje následující funkce:
  

	A = M.Matrix([[5,5],[5,5]] - Načtení pomocí 2D seznamu
	 	
	A = M.Matrix.from_list([1,2,3,4,5,6,7,8,9],3,3)   - vytvoří matici s hodnotami [[1,2,3],[4,5,6],[7,8,9]] - Načtení pomocí 1D seznamu s údajem kolik má matice mít řádků a sloupc
	
	A.type() -  vypíše informace jestli je matice: nulová, jednotková, čtvercová, diagonální, symetrická, regulární
	
	A.nulova()  - vrací True nebo False pokud je matice nulová nebo ne
 
	A.jednotkova()  - vrací True nebo False pokud je matice jednotková nebo ne
 
	A.ctvercova()  - vrací True nebo False pokud je matice čtvercová nebo ne
 
	A.diagonálni()  - vrací True nebo False pokud je matice diagonální nebo ne
 
	A.symetricka()  - vrací True nebo False pokud je matice symetrická nebo ne
 
	A.regularni()  - vrací True nebo False pokud je matice regulární nebo ne
 
	A.ortogonalni() – vrací True nebo False pokud je matice ortogonální nebo ne
 
	A.pozitivne_definitni() – vrací True nebo False pokud je matice pozitivně definitní nebo ne
 
	A.transpozice() – vrací matici transponovanou k A
 
	A + B vrací součet matic A a B (lze použít také A.__add__(B))
 
	A – B vrací odečtení matic A a B(lze použít také A.__sub__(B))
 
	A*B vrací maticový součin matic A a B (lze použít také A.__mul__(B))
 
	A.kron(B) - vrací kroneckerovský součin matic A a B
 
	A.Gauss() – vrací matici A po provedení Gaussovy eliminace (destruktivní pro matici A)
 
	A.Gauss_Jordan() vrací matici A po provedení Gauss-Jordanovy eliminaci (destruktivní pro matici A)  
 
	A.hodnost() – vrací hodnost matice A
 
	A.dimenze_jadra() – vrací dimenzi jádra matice A
 
	A.reseni(prava_strana) – vrací vektor x jakožto řešení rovnice A*x = prava_strana, kde pravá strana je nějaký vektor
 
	A-Ker() – vrací jádro matice A
 
	A.determinant() – vrací determinant matice A
 
	A.inverz() – vrací matici inverzní k A
 
	v1.skalarni_soucin(v2) – vrací skalární součin 2 vektorů (Maticí s 1 sloupcem)
 
	v1.norma() – vrací normu vektoru v1
 
	A.QR() – vrací 2 matice Q a R z QR rozkladu matice A (použití: Q, R = A.QR())
 
	A.vlastni_cisla() - vrací reálná vlastní čísla matice A
 	
	A.vlastni_vektory() – vrací reálné vlastní vektory matice A
 
	A.Choleskeho_rozklad() – vrací matici L z cholevskeho rozkladu matice A


Testování:
-	Pro knihovnu jsem vytvořil testovací kód na jednotkové testy pomocí knihovny unittest
-	V tomto kódu jsem vytvořil několik matic a vektorů s různými vlastnostmi pro otestování všech funkcí
-	Importoval jsem také knihovnu numpy, díky které mohu zároveň testovat i validovat moji knihovnu
-	Některé funkce, které jsem vytvořil pro můj class Matrix v numpy nejsou(nebo jsem je alespoň nedokázal najít – např Gaussova eliminace), ale díky matematickému spojení jednotlivých funkci jsem otestoval všechny funkce, že fungují správně. Tedy kupříkladu zmiňovaná Gaussova eliminace je testována třeba při testu řešení rovnice A*x = b 
-	Všechny funkce tedy fungují správně a jsou validovány s funkcemi v numpy

 

Doplňující informace:
-	Do knihovny jsem se rozhodl přidat práci s komplexními čísly, takže knihovna dokáže pracovat s komplexními čísly i komplexními maticemi, kromě výpočtu vlastních čísel a vektorů
-	Snažil jsem se knihovnu udělat uživatelsky pohodlnou, proto jsem zvolil zaokrouhlování výsledků po pythonských chybách a zaokrouhlení výsledků na 3 desetinná místa(např u determinantu). U ukládání hodnot matic ve výpočtech se všechny hodnoty ukládají s maximální přesností, ale při výpisu a reprezentaci matice se vypíšou opět zaokrouhlené hodnoty na 3 desetinná místa, pro lepší přehlednost. Výpis není destruktivní pro původní přesné hodnoty a uloženy stále zůstávají s maximální přesností pro další počítání.
-	Samozřejmostí pro všechny funkce je ošetřený vstup – pokud uživatel zadá neplatný vstup = chce násobit dvě nekompatibilní matice, provádí skalární součin matice která není vektor a podobně, kód vyhodí error s informací co je na vstupu špatně
