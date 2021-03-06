package com.Main;

import java.util.ArrayList;
import Classes.Mesh;
import Tools.Matrix;
import Tools.Vector;
import java.util.ArrayList;
import static Enums.Sizes.*;
import static Tools.Math_Tools.*;
import static Tools.Sel.*;
import static Tools.Tools.*;
public class Main {

    public static void main(String[] args) {
        System.out.println("POLYGLOT SIDEQUEST");
        //Se guarda el nombre del archivo obtenido de los argumentos pasados al programa
       String filename = args[0];
        //Se instancian los vectores y matrices a utilizar
        ArrayList<Matrix> localKs = new ArrayList<Matrix>();
        ArrayList<Vector> localbs = new ArrayList<Vector>();
        Matrix K = new Matrix();
        Vector b = new Vector();
        Vector T = new Vector();

        //Presentacion al programa
        System.out.println("IMPLEMENTACION DEL METODO DE LOS ELEMENTOS FINITOS EN 3D");
        System.out.println("\t- POLIGLOT SIDEQUEST\n"+"\t- 3 DIMENSIONES");
        System.out.println("\t- FUNCIONES DE FORMA LINEALES\n"+"\t- PESOS DE GALERKIN");
        System.out.println("\t- MALLA CUADRATICA DE TETAEDRO IRREGULAR");
        System.out.println("*************************************************************************");

        //Se instancia el objeto de la clase Mesh
        Mesh m = new Mesh();
        //Se leen los datos del archivo .dat
        leerMallayCondiciones(m, filename);
        System.out.println("Datos obtenidos correctamente\n***********************");

        //Se crean los sitemas locales y se muestran
        crearSistemasLocales(m, localKs, localbs);
        //showKs(localKs); showbs(localbs);
        System.out.println("*******************************");

        //Las matrices K y b se llenan de 0 y posteriormente se realiza el ensamblaje
        zeroes(K, m.getSize(NODES.ordinal()));
        zeroes(b, m.getSize(NODES.ordinal()));
        ensamblaje(m, localKs, localbs, K, b);

        System.out.println("*******************************");

        //Se aplica la condicion de Neumann
        applyNeumann(m, b);
        System.out.println("*******************************");

        //Se aplica la condicion de dirichlet
        applyDirichlet(m,K,b);
        System.out.println("*******************************");

        //Se calcula el resultado del SEL
        zeroes(T, b.size());
        calculate(K,b,T);

        //Se crea un archivo con los resultados.
        writeResults(m,T,filename);

    }
}
