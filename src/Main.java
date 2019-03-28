public class Main {

    public static void main(String[] args) {

        // https://planetcalc.com/3324/ za test
        double[][] mat = new double[][]{{6, 24, 1}, {13, 16, 10}, {20, 17, 15}};
        int[][] modInv = modularInverse(mat, 26);

        for (int i = 0; i < modInv.length; i++) {

            for (int j = 0; j < modInv.length; j++) {

                System.out.print(modInv[i][j] + " ");
            }

            System.out.println();
        }
    }

    public static int[][] modularInverse(double[][] mat, int mod) {

        // 1. determinanta
        int det = modulus((int) determinant(mat), mod);

        // TODO prej preveri če je determinanta prava ( != 0 in je tuje z mod)
        // 2. zdaj je treba zracunati enacbo (DET * x - 1) mod (mod) == 0, dobiti moramo x
        int x = getX(det, mod);

        // 3. adjungirana oblika
        double[][] adj = adjugateForm(mat);

        // 4. izracun mod inv matrike z vsemi podatki
        return calcModInv(x, adj, mod);
    }

    // zracuna koncnen modularni invetz iz vseh podatkov ki jih rabimo: adjungiana oblika, x in modulus
    // modularni inverz dobimo tako da vpomnozimo x z adjungirano matriko, nato pa element se modulamo z mod-om
    public static int[][] calcModInv(int x, double[][] adj, int mod) {

        int[][] modInv = new int[adj.length][adj.length];

        for (int i = 0; i < adj.length; i++) {

            for (int j = 0; j < adj.length; j++) {

                modInv[i][j] = modulus((int) (x * adj[i][j]), mod);
            }
        }

        return modInv;
    }

    // zracuna n mod (mod)
    public static int modulus(int n, int mod) {

        if (n >= 0) {
            return n % mod;
        } else {
            return -(Math.abs(n) % mod) + mod;
        }
    }

    // izracuna enacbo (N * x - 1) mod (mod) == 0 po bruteforce metodi
    public static int getX(int N, int mod) {

        int x = 0;

        while (true) {

            if (modulus(N * x - 1, mod) == 0) {

                return x;
            }
            x++;
        }
    }

    // adjungirana matrika matrike mat zgleda tako, da:
    // prekrijemo i-to vrstico in -jti stolpevc, od tega kar ostane zracunamo determinanto
    // determinanto vstavimo na i, j -to mesto in dodamo rpedznak, ki ima boliko sahovnice, zacne se ze +
    // dano matriko je treba se transponirat
    public static double[][] adjugateForm(double[][] mat) {

        // adjungirana od enojne matrike je kar ta matrika
        if (mat.length == 1) {
            return mat;
        }

        double[][] adjugateMatrix = new double[mat.length][mat.length];
        double[][] tmpMatrix;

        for (int i = 0; i < mat.length; i++) {

            for (int j = 0; j < mat.length; j++) {

                // prekrij i-to vrstico in j-ti stolpec in naredi matriko iz ostalih elementov
                tmpMatrix = matrixWithoutRowAndColumn(i, j, mat);

                // pridobi determinanto te matrike in jo shrani na pravo mesto
                // dodaj -, če potrebno: v obliki šahovnice
                /*
                + - + - + - ..
                - + - + ..
                + - ..
                 */
                // j, i je zato, ker se mora matrika transponirat
                adjugateMatrix[j][i] = ((i + j) % 2 == 0 ? 1 : -1) * determinant(tmpMatrix);
            }
        }

        return adjugateMatrix;
    }

    // zračunaj amtriko brez i-te vrstice in j-tega stolpca
    public static double[][] matrixWithoutRowAndColumn(int row, int col, double[][] mat) {

        double[][] tmpMatrix = new double[mat.length - 1][mat.length - 1];
        int tmpCol = 0, tmpRow = 0;

        for (int i = 0; i < mat.length; i++) {

            for (int j = 0; j < mat.length; j++) {

                if (!(i == row || j == col)) {

                    tmpMatrix[tmpRow][tmpCol] = mat[i][j];
                    tmpCol++;
                }
            }

            tmpCol = 0;

            if (i != row) {

                tmpRow++;
            }
        }

        return tmpMatrix;
    }

    // zracunaj determinanto matirke
    public static double determinant(double[][] mat) {

        // za enojne matrike je kar vrednost v matriki
        if (mat.length == 1) {

            return mat[0][0];
        }

        // za dvojne matrike je enacba ad - bc
        if (mat.length == 2) {

            return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
        }

        // prekrita matrika
        double[][] tmp = new double[mat.length - 1][mat.length - 1];
        double det = 0;

        // za matrike nad 2 pa:
        // - gremo cez vsaki element na prvi vrstici
        for (int i = 0; i < mat.length; i++) {

            // - pri vsakem elementu na i,j ti poziciji pokrijemo vrstico 0 in stolpec j
            // in naredimo matriko iz preostalih elementov (ima eno vrstico in stolpec manj)
            tmp = matrixWithoutRowAndColumn(i, 0, mat);

            // zdaj pa izračunamo determinanto te podmatrike rekurzivno in pomnozimo s ternutnim elementon
            // ter s predznakom ki se menjuje, zacne se z +
            // vse skupaj dodamo v koncno determinanto
            det += mat[i][0] * (i % 2 == 0 ? 1 : -1) * determinant(tmp);
        }

        return det;
    }
}
