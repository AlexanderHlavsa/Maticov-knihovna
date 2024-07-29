import unittest
import Maticova_knihovna as M
import numpy as np

A = M.Matrix.from_list([1,1,0,2,-1,1,0,1,1,1,0,3],4,3)
A1 = np.array([[1, 1, 0], [2, -1, 1], [0, 1, 1], [1, 0, 3]])

B = M.Matrix.from_list([9,4,7,6,5,4,3,2,1],3,3)
B1= np.array([[9, 4, 7], [6, 5, 4], [3, 2, 1]])

C = M.Matrix.from_list([2,2,1,1,1,1],3,2)
C1= np.array([[2, 2], [1, 1], [1, 1]])

D = M.Matrix.from_list([1,6,1,2,8,9,2,3,2],3,3)
D1 = np.array([[1, 6, 1], [2, 8, 9], [2, 3, 2]])

J = M.Matrix.from_list([1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],4,4)
J1 = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])

N = M.Matrix.from_list([0,0,0,0,0,0,0,0,0],3,3)
N1 = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])

P = M.Matrix.from_list([1,2,3,-1,2,1,0,0,1,2,0,0,2,1,2,1,0,1,0,0,0],3,7)
P1= np.array([[1, 2, 3, -1, 2, 1, 0], [0, 1, 2, 0, 0, 2, 1], [2, 1, 0, 1, 0, 0, 0]])

S = M.Matrix([[1,2,3,4],[2,1,0,0],[3,0,1,5],[4,0,5,1]])
S1 = np.array([[1,2,3,4],[2,1,0,0],[3,0,1,5],[4,0,5,1]])

G = M.Matrix.from_list([5,5,5,5,5,5,5,5,5,7,5,5,5,5,6],5,3)
G1= np.array([[5, 5, 5], [5, 5, 5], [5, 5, 5], [7, 5, 5], [5, 5, 6]])

H = M.Matrix([[4,-2,4],[-2,10,1],[4,1,6]])
H1 = np.array([[4,-2,4],[-2,10,1],[4,1,6]])

I = M.Matrix([[1,-3,2],[-3,13,2],[2,2,21]])
I1 = np.array([[1,-3,2],[-3,13,2],[2,2,21]])

L = M.Matrix([[5,5],[5,25]])
L1 = np.array([[5,5],[5,25]])

K = M.Matrix([[5 + 1j, 6, 1],[2, 8 + 2j, 9],[2, 3, 2 - 1j]])
K1 = np.array([[5 + 1j, 6, 1],[2, 8 + 2j, 9],[2, 3, 2 - 1j]])

prava_strana1 = M.Matrix.from_list([2,0,5],3,1)
prava_strana1_1 = np.array([[2], [0], [5]])

prava_strana2 = M.Matrix([[6],[-3],[5],[1]])
prava_strana2_1 = np.array([[6],[-3],[5],[1]])

prava_strana3 = M.Matrix([[-10],[-10],[-10],[-12],[-10]])
prava_strana3_1 = np.array([[-10],[-10],[-10],[-12],[-10]])

#testovani
class TestMatrix(unittest.TestCase):
    #transpozice
    def test_transpoziceA(self):
        self.assertTrue((A.transpozice().data == A1.transpose()).all())
    def test_transpoziceB(self):
        self.assertTrue((B.transpozice().data == B1.transpose()).all())
    def test_transpoziceC(self):
        self.assertTrue((C.transpozice().data == C1.transpose()).all())
    def test_transpoziceD(self):
        self.assertTrue((D.transpozice().data == D1.transpose()).all())
    def test_transpoziceN(self):
        self.assertTrue((N.transpozice().data == N1.transpose()).all())
    def test_transpoziceP(self):
        self.assertTrue((P.transpozice().data == P1.transpose()).all())
    def test_transpoziceG(self):
        self.assertTrue((G.transpozice().data == G1.transpose()).all())
    def test_transpoziceS(self):
        self.assertTrue((S.transpozice().data == S1.transpose()).all())
    def test_transpoziceJ(self):
        self.assertTrue((J.transpozice().data == J1.transpose()).all())
    def test_transpoziceK(self):
        self.assertTrue((K.transpozice().data == K1.transpose()).all())
    #type
    def test_typeA(self):
        self.assertEqual(A.type(),"")
    def test_typeB(self):
        self.assertEqual(B.type(),"ctvercova, regularni")
    def test_typeC(self):
        self.assertEqual(C.type(),"")
    def test_typeD(self):
        self.assertEqual(D.type(),"ctvercova, regularni")
    def test_typeJ(self):
        self.assertEqual(J.type(),"jednotkova, diagonalni, symetricka, ctvercova, regularni")
    def test_typeS(self):
        self.assertEqual(S.type(),"symetricka, ctvercova, regularni")
    def test_typeN(self):
        self.assertEqual(N.type(),"nulova, ctvercova")
    def test_typeG(self):
        self.assertEqual(G.type(),"")
    def test_typeP(self):
        self.assertEqual(P.type(),"")

    #determinant, pouze u ctvercovych
    def test_determinantB(self):
        self.assertAlmostEqual(B.determinant(),np.linalg.det(B1))
    def test_determinantD(self):
        self.assertAlmostEqual(D.determinant(),np.linalg.det(D1))
    def test_determinantJ(self):
        self.assertAlmostEqual(J.determinant(),np.linalg.det(J1))
    def test_determinantN(self):
        self.assertAlmostEqual(N.determinant(),np.linalg.det(N1))
    def test_determinantS(self):
        self.assertAlmostEqual(S.determinant(),np.linalg.det(S1))
    def test_determinantK(self):
        self.assertAlmostEqual(K.determinant(),np.linalg.det(K1))
    #soucet
    def test_soucetB_D(self):
        self.assertTrue(((B+D).data == (D+B).data == B1+D1).all())
    def test_soucetB_N(self):
        self.assertTrue(((B+N).data == (N+B).data == B1+N1).all())
    def test_soucetN_D(self):
        self.assertTrue(((N+D).data == (D+N).data == N1+D1).all())
    def test_soucetJ_S(self):
        self.assertTrue(((J+S).data == (S+J).data == J1+S1).all())
    def test_soucetB_K(self):
        self.assertTrue(((B+K).data == (K+B).data == B1+K1).all())
    #odecitani
    def test_odecitaniB_D(self):
        self.assertTrue(((B-D).data == (-(D-B)).data == B1-D1).all())
    def test_odecitaniB_N(self):
        self.assertTrue(((B-N).data == (-(N-B)).data == B1-N1).all())
    def test_odecitaniN_D(self):
        self.assertTrue(((N-D).data == (-(D-N)).data == N1-D1).all())
    def test_odecitaniJ_S(self):
        self.assertTrue(((J-S).data == (-(S-J)).data == J1-S1).all())

    #nasobeni
    def test_soucinJ_A(self):
        self.assertTrue(((J*A).data == J1.dot(A1)).all())
    def test_soucinA_B(self):
        self.assertTrue(((A*B).data == A1.dot(B1)).all())
    def test_soucinA_C(self):
        self.assertTrue(((A*C).data == A1.dot(C1)).all())
    def test_soucinB_C(self):
        self.assertTrue(((B*C).data == B1.dot(C1)).all())
    def test_soucinG_P(self):
        self.assertTrue(((G*P).data == G1.dot(P1)).all())
    def test_soucinG_N(self):
        self.assertTrue(((G*N).data == G1.dot(N1)).all())
    def test_soucinD_P(self):
        self.assertTrue(((D*P).data == D1.dot(P1)).all())
    def test_soucinD_K(self):
        self.assertTrue(((D*K).data == D1.dot(K1)).all())
    def test_soucinB_prava_strana(self):
        self.assertTrue(((B*prava_strana1).data == B1.dot(prava_strana1_1)).all)

    #kronecker
    def test_kronA_B(self):
        self.assertTrue((A.kron(B).data == np.kron(A1,B1)).all())
    def test_kronB_C(self):
        self.assertTrue((B.kron(C).data == np.kron(B1,C1)).all())
    def test_kronC_D(self):
        self.assertTrue((C.kron(D).data == np.kron(C1,D1)).all())
    def test_kronD_J(self):
        self.assertTrue((D.kron(J).data == np.kron(D1,J1)).all())
    def test_kronJ_N(self):
        self.assertTrue((J.kron(N).data == np.kron(J1,N1)).all())
    def test_kronN_P(self):
        self.assertTrue((N.kron(P).data == np.kron(N1,P1)).all())
    def test_kronP_S(self):
        self.assertTrue((P.kron(S).data == np.kron(P1,S1)).all())
    def test_kronS_G(self):
        self.assertTrue((S.kron(G).data == np.kron(S1,G1)).all())
    def test_kronS_K(self):
        self.assertTrue((S.kron(K).data == np.kron(S1,K1)).all())
    def test_kronG_prav(self):
        self.assertTrue((G.kron(prava_strana1).data == np.kron(G1,prava_strana1_1)).all())
    
    #hodnost
    def test_hodnostA(self):
        self.assertEqual(A.hodnost(), np.linalg.matrix_rank(A1))
    def test_hodnostB(self):
        self.assertEqual(B.hodnost(), np.linalg.matrix_rank(B1))
    def test_hodnostC(self):
        self.assertEqual(C.hodnost(), np.linalg.matrix_rank(C1))
    def test_hodnostD(self):
        self.assertEqual(D.hodnost(), np.linalg.matrix_rank(D1))
    def test_hodnostJ(self):
        self.assertEqual(J.hodnost(), np.linalg.matrix_rank(J1))
    def test_hodnostN(self):
        self.assertEqual(N.hodnost(), np.linalg.matrix_rank(N1))
    def test_hodnostP(self):
        self.assertEqual(P.hodnost(), np.linalg.matrix_rank(P1))
    def test_hodnostS(self):
        self.assertEqual(S.hodnost(), np.linalg.matrix_rank(S1))
    def test_hodnostG(self):
        self.assertEqual(G.hodnost(), np.linalg.matrix_rank(G1))
    def test_hodnostK(self):
        self.assertEqual(K.hodnost(), np.linalg.matrix_rank(K1))
    def test_hodnostPrava_str(self):
        self.assertEqual(prava_strana1.hodnost(), np.linalg.matrix_rank(prava_strana1_1))

    #reseni s pravou stranou = zaroven i kontrola funkci Gauss,Gauss_Jordan a Jordan
    def test_reseniB_prava_strana1(self):
        X = np.array(B.reseni(prava_strana1).data)
        self.assertTrue(np.allclose(np.dot(B1,X),prava_strana1_1))
    def test_reseniC_prava_strana1(self):
        with self.assertRaises(ValueError):
            C.reseni(prava_strana1)
    def test_reseniD_prava_strana1(self):
        X = np.array(D.reseni(prava_strana1).data)
        self.assertTrue(np.allclose(np.dot(D1,X),prava_strana1_1))
    def test_reseniN_prava_strana1(self):
        with self.assertRaises(ValueError):
            N.reseni(prava_strana1)
    def test_reseniP_prava_strana1(self):
        X = np.array(P.reseni(prava_strana1).data)
        self.assertTrue(np.allclose(np.dot(P1,X),prava_strana1_1))
    def test_reseniA_prava_strana2(self):
        X = np.array(A.reseni(prava_strana2).data)
        self.assertTrue(np.allclose(np.dot(A1,X),prava_strana2_1))
    def test_reseniJ_prava_strana2(self):
        X = np.array(J.reseni(prava_strana2).data)
        self.assertTrue(np.allclose(np.dot(J1,X),prava_strana2_1))
    def test_reseniS_prava_strana2(self):
        X = np.array(S.reseni(prava_strana2).data)
        self.assertTrue(np.allclose(np.dot(S1,X),prava_strana2_1))
    def test_reseniG_prava_strana3(self):
        X = np.array(G.reseni(prava_strana3).data)
        self.assertTrue(np.allclose(np.dot(G1,X),prava_strana3_1))

    #Jadro
    def test_jadroA(self):
        jadro = A.Ker()
        hodnost = A.hodnost()
        if A.sloupce != hodnost:           
            self.assertEqual((A.sloupce)-hodnost,len(jadro))
            nul = np.array([[0]for _ in range(A.radky)])
            for vektor in jadro:
                self.assertTrue(np.allclose(np.dot(A1,np.array(vektor.data)),nul))
    def test_jadroB(self):
        jadro = B.Ker()
        hodnost = B.hodnost()
        if B.sloupce != hodnost:           
            self.assertEqual((B.sloupce)-hodnost,len(jadro))
            nul = np.array([[0]for _ in range(B.radky)])
            for vektor in jadro:
                self.assertTrue(np.allclose(np.dot(B1,np.array(vektor.data)),nul))
    def test_jadroC(self):
        jadro = C.Ker()
        hodnost = C.hodnost()
        if C.sloupce != hodnost:           
            self.assertEqual((C.sloupce)-hodnost,len(jadro))
            nul = np.array([[0]for _ in range(C.radky)])
            for vektor in jadro:
                self.assertTrue(np.allclose(np.dot(C1,np.array(vektor.data)),nul))
    def test_jadroD(self):
        jadro = D.Ker()
        hodnost = D.hodnost()
        if D.sloupce != hodnost:           
            self.assertEqual((D.sloupce)-hodnost,len(jadro))
            nul = np.array([[0]for _ in range(D.radky)])
            for vektor in jadro:
                self.assertTrue(np.allclose(np.dot(D1,np.array(vektor.data)),nul))
    def test_jadroJ(self):
        jadro = J.Ker()
        hodnost = J.hodnost()
        if J.sloupce != hodnost:           
            self.assertEqual((J.sloupce)-hodnost,len(jadro))
            nul = np.array([[0]for _ in range(J.radky)])
            for vektor in jadro:
                self.assertTrue(np.allclose(np.dot(J1,np.array(vektor.data)),nul))
    def test_jadroN(self):
        jadro = N.Ker()
        hodnost = N.hodnost()
        if N.sloupce != hodnost:           
            self.assertEqual((N.sloupce)-hodnost,len(jadro))
            nul = np.array([[0]for _ in range(N.radky)])
            for vektor in jadro:
                self.assertTrue(np.allclose(np.dot(N1,np.array(vektor.data)),nul))
    def test_jadroP(self):
        jadro = P.Ker()
        hodnost = P.hodnost()
        if P.sloupce != hodnost:           
            self.assertEqual((P.sloupce)-hodnost,len(jadro))
            nul = np.array([[0]for _ in range(P.radky)])
            for vektor in jadro:
                self.assertTrue(np.allclose(np.dot(P1,np.array(vektor.data)),nul))
    def test_jadroS(self):
        jadro = S.Ker()
        hodnost = S.hodnost()
        if S.sloupce != hodnost:           
            self.assertEqual((S.sloupce)-hodnost,len(jadro))
            nul = np.array([[0]for _ in range(S.radky)])
            for vektor in jadro:
                self.assertTrue(np.allclose(np.dot(S1,np.array(vektor.data)),nul))
    def test_jadroG(self):
        jadro = G.Ker()
        hodnost = G.hodnost()
        if G.sloupce != hodnost:           
            self.assertEqual((G.sloupce)-hodnost,len(jadro))
            nul = np.array([[0]for _ in range(G.radky)])
            for vektor in jadro:
                self.assertTrue(np.allclose(np.dot(G1,np.array(vektor.data)),nul))
    def test_jadroK(self):
        jadro = K.Ker()
        hodnost = K.hodnost()
        if K.sloupce != hodnost:           
            self.assertEqual((K.sloupce)-hodnost,len(jadro))
            nul = np.array([[0]for _ in range(K.radky)])
            for vektor in jadro:
                self.assertTrue(np.allclose(np.dot(K1,np.array(vektor.data)),nul))

    #inverzni matice
    def test_inverzA(self):
        with self.assertRaises(ValueError):
            A.inverz()
    def test_inverzC(self):
        with self.assertRaises(ValueError):
            C.inverz()
    def test_inverzP(self):
        with self.assertRaises(ValueError):
            P.inverz()
    def test_inverzG(self):
        with self.assertRaises(ValueError):
            G.inverz()
    def test_inverzB(self):
        self.assertTrue(np.allclose(B.inverz().data, np.linalg.inv(B1)))
    def test_inverzD(self):
        self.assertTrue(np.allclose(D.inverz().data, np.linalg.inv(D1)))
    def test_inverzJ(self):
        self.assertTrue(np.allclose(J.inverz().data, np.linalg.inv(J1)))
    def test_inverzS(self):
        self.assertTrue(np.allclose(S.inverz().data, np.linalg.inv(S1)))
    def test_inverzK(self):
        self.assertTrue(np.allclose(K.inverz().data, np.linalg.inv(K1)))

    #QR rozklad
    def test_QR_A(self):
        Q,R = A.QR()
        self.assertTrue(Q.ortogonalni())
        self.assertTrue(R.horni_trojuhelnikova())
        self.assertTrue(np.allclose((Q*R).data,A1))
    def test_QR_B(self):
        Q,R = B.QR()
        self.assertTrue(Q.ortogonalni())
        self.assertTrue(R.horni_trojuhelnikova())
        self.assertTrue(np.allclose((Q*R).data,B1))
    def test_QR_D(self):
        Q,R = D.QR()
        self.assertTrue(Q.ortogonalni())
        self.assertTrue(R.horni_trojuhelnikova())
        self.assertTrue(np.allclose((Q*R).data,D1))
    def test_QR_J(self):
        Q,R = J.QR()
        self.assertTrue(Q.ortogonalni())
        self.assertTrue(R.horni_trojuhelnikova())
        self.assertTrue(np.allclose((Q*R).data,J1))
    def test_QR_S(self):
        Q,R = S.QR()
        self.assertTrue(Q.ortogonalni())
        self.assertTrue(R.horni_trojuhelnikova())
        self.assertTrue(np.allclose((Q*R).data,S1))
    def test_QR_G(self):
        Q,R = G.QR()
        self.assertTrue(Q.ortogonalni())
        self.assertTrue(R.horni_trojuhelnikova())
        self.assertTrue(np.allclose((Q*R).data,G1))
    def test_QR_K(self):
        Q,R = K.QR()
        self.assertTrue(Q.ortogonalni())
        self.assertTrue(R.horni_trojuhelnikova())
        self.assertTrue(np.allclose((Q*R).data,K1))
    def test_QR_C(self):
        with self.assertRaises(ValueError):
            C.QR()
    def test_QR_N(self):
        with self.assertRaises(ValueError):
            N.QR()
    def test_QR_P(self):
        with self.assertRaises(ValueError):
            P.QR()

    #vlastni cisla a vektory
    def test_vlastni_c_B(self):
        for v_cislo in B.vlastni_cisla():
            self.assertTrue(any(np.allclose(v_cislo,vlastni_cislo) for vlastni_cislo in np.linalg.eig(B1)[0]))
    def test_vlastni_c_D(self):
        for v_cislo in D.vlastni_cisla():
            self.assertTrue(any(np.allclose(v_cislo,vlastni_cislo) for vlastni_cislo in np.linalg.eig(D1)[0]))
    def test_vlastni_c_J(self):
        for v_cislo in J.vlastni_cisla():
            self.assertTrue(any(np.allclose(v_cislo,vlastni_cislo) for vlastni_cislo in np.linalg.eig(J1)[0]))
    def test_vlastni_c_N(self):
        for v_cislo in N.vlastni_cisla():
            self.assertTrue(any(np.allclose(v_cislo,vlastni_cislo) for vlastni_cislo in np.linalg.eig(N1)[0]))
    def test_vlastni_c_S(self):
        for v_cislo in S.vlastni_cisla():
            self.assertTrue(any(np.allclose(v_cislo,vlastni_cislo) for vlastni_cislo in np.linalg.eig(S1)[0]))
    def test_vlastni_c_K(self):
        for v_cislo in K.vlastni_cisla():
            self.assertTrue(any(np.allclose(v_cislo,vlastni_cislo) for vlastni_cislo in np.linalg.eig(K1)[0]))
    def test_vlastni_c_A(self):
        with self.assertRaises(ValueError):
            A.vlastni_cisla()
    def test_vlastni_c_C(self):
        with self.assertRaises(ValueError):
            C.vlastni_cisla()
    def test_vlastni_c_P(self):
        with self.assertRaises(ValueError):
            P.vlastni_cisla()
    def test_vlastni_c_G(self):
        with self.assertRaises(ValueError):
            G.vlastni_cisla()

    #vlastni_vektory se testuji primo v kodu
    
    def test_cholesky_H(self):
        self.assertTrue(np.allclose(H.Choleskeho_rozklad().data,np.linalg.cholesky(H1)))
    def test_cholesky_I(self):
        self.assertTrue(np.allclose(I.Choleskeho_rozklad().data,np.linalg.cholesky(I1)))
    def test_cholesky_L(self):
        self.assertTrue(np.allclose(L.Choleskeho_rozklad().data,np.linalg.cholesky(L1)))
if __name__ == '__main__':
    unittest.main()
