package Tools;
import Classes.*;
import Enums.Parameters;
import Enums.Sizes;

import java.util.ArrayList;
import static Tools.Math_Tools.*;
public class Sel {

    private Sel(){};

    //Muestra las Ks
    public static void showKs(ArrayList<Matrix> Ks){
        for (int i = 0; i < Ks.size() ; i++) {
            System.out.println("K del elemento"+(i+1));
            Ks.get(i).Show();
            System.out.println("**********************************");
        }
    }
    //obtencion del valor C1

    public static float CreateC1 (int ind, Mesh m){
            Element el = m.getElement(ind);
            el= m.getElement(ind);
            Node n1 = m.getNode(el.getNode1()-1);
            Node n2 = m.getNode(el.getNode2()-1);



            return (float) ((1f)/Math.pow((n2.getX()-n1.getX()),2));


    }

    //obtencion del C2
    public static float CreateC2 (int ind, Mesh m){
        Element el = new Element();
        el= m.getElement(ind);
        Node n1 = m.getNode(el.getNode1()-1);
        Node n2 = m.getNode(el.getNode2()-1);
        Node n8 = m.getNode(el.getNode8()-1);


        return (float) ((1f)/(n2.getX() - n1.getX()))*((4*n1.getX())+(4*n2.getX())-(8*n8.getX()));


    }
    //crenado matriz μ
    public static void calculateLocalμ(float c1, float c2, Matrix μ){

        float A,B,C,D,E,F,G,H,I,J,K;


        A= (float) ((float) ((float) ((float) (-(1f)/(192*Math.pow(c2,2)))*Math.pow(4*c1-c2,4))-((1f)/(24*c2))*Math.pow(4*c1-c2,3))-
                ((1f)/(3840*Math.pow(c2,3)))*Math.pow(4*c1-c2,5)+((1f)/(3840*Math.pow(c2,3)))*Math.pow(4*c1-3*c2,5));

        B= (float) ((float) ((float) ((float) (-(1f)/(192*Math.pow(c2,2)))*Math.pow(4*c1+c2,4))+((1f)/(24*c2))*Math.pow(4*c1+c2,3))+
                ((1f)/(3840*Math.pow(c2,3)))*Math.pow(4*c1+c2,5)-((1f)/(3840*Math.pow(c2,3)))*Math.pow(4*c1-3*c2,5));

        C= (float) ((4/15)*Math.pow(c2,2));

        D= (float) ((float) ((float) ((float) ((float) ((float) ((float) ((float) ((float) ((1f)/(192*Math.pow(c2,2)))*Math.pow(4*c2-c1,4))-((1f)/(3840*Math.pow(c2,3)))*Math.pow(4*c2-c1,5))+
                                                        ((1f)/(7680*Math.pow(c2,3)))*Math.pow(4*c2+8*c1,5)-((7f)/(7680*Math.pow(c2,3)))*Math.pow(4*c2-8*c1,5))+ ((1f)/(768*Math.pow(c2,3)))*Math.pow(-8*c1,5))-
                                                (c1)/(96*Math.pow(c2,3)))*Math.pow((4*c2)-8*c1,4))+((2*c1-1)/192*Math.pow(c2,3)))*Math.pow(-8*c1,4));

        E = (float) ((float) ((8/3)*Math.pow(c1,2))+ (1/30)*Math.pow(c2,2));

        F= (float) (((2/3)*(c1*c2)) - (1/30)*Math.pow(c2,2));

        G = (float) ((-16/3)*(Math.pow(c1,2)) - (4/3)*(c1*c2) - (2/15)*Math.pow(c2,2)) ;

        H = (float) (((2/3)*(c1*c2)) + (1/30)*Math.pow(c2,2));

        I =  (float) ((float) (-(16/3)*Math.pow(c1,2))-(2/3)*Math.pow(c2,2));

        J = (float) (- (2/15)*Math.pow(c2,2)) ;

        K =  -(4/3)*(c1*c2);
        //seteando valores de μ
        μ.get(0).set(0, A);
        μ.get(0).set(1, E);
        μ.get(0).set(2, 0f);
        μ.get(0).set(3, 0f);
        μ.get(0).set(4, -F);
        μ.get(0).set(5, 0f);
        μ.get(0).set(6, -F);
        μ.get(0).set(7, G);
        μ.get(0).set(8, F);
        μ.get(0).set(9, F);

        μ.get(1).set(0, E);
        μ.get(1).set(1, B);
        μ.get(1).set(2, 0f);
        μ.get(1).set(3, 0f);
        μ.get(1).set(4, -H);
        μ.get(1).set(5, 0f);
        μ.get(1).set(6, -H);
        μ.get(1).set(7, I);
        μ.get(1).set(8, H);
        μ.get(1).set(9, H);

        μ.get(2).set(0, 0f);
        μ.get(2).set(1, 0f);
        μ.get(2).set(2, 0f);
        μ.get(2).set(3, 0f);
        μ.get(2).set(4, 0f);
        μ.get(2).set(5, 0f);
        μ.get(2).set(6, 0f);
        μ.get(2).set(7, 0f);
        μ.get(2).set(8, 0f);
        μ.get(2).set(9, 0f);

        μ.get(3).set(0, 0f);
        μ.get(3).set(1, 0f);
        μ.get(3).set(2, 0f);
        μ.get(3).set(3, 0f);
        μ.get(3).set(4, 0f);
        μ.get(3).set(5, 0f);
        μ.get(3).set(6, 0f);
        μ.get(3).set(7, 0f);
        μ.get(3).set(8, 0f);
        μ.get(3).set(9, 0f);

        μ.get(4).set(0, -F);
        μ.get(4).set(1, -H);
        μ.get(4).set(2, 0f);
        μ.get(4).set(3, 0f);
        μ.get(4).set(4, C);
        μ.get(4).set(5, 0f);
        μ.get(4).set(6, J);
        μ.get(4).set(7, -K);
        μ.get(4).set(8, -C);
        μ.get(4).set(9, -J);

        μ.get(5).set(0, 0f);
        μ.get(5).set(1, 0f);
        μ.get(5).set(2, 0f);
        μ.get(5).set(3, 0f);
        μ.get(5).set(4, 0f);
        μ.get(5).set(5, 0f);
        μ.get(5).set(6, 0f);
        μ.get(5).set(7, 0f);
        μ.get(5).set(8, 0f);
        μ.get(5).set(9, 0f);

        μ.get(6).set(0, -F);
        μ.get(6).set(1, -H);
        μ.get(6).set(2, 0f);
        μ.get(6).set(3, 0f);
        μ.get(6).set(4, J);
        μ.get(6).set(5, 0f);
        μ.get(6).set(6, C);
        μ.get(6).set(7, -K);
        μ.get(6).set(8, -J);
        μ.get(6).set(9, -C);

        μ.get(7).set(0, G);
        μ.get(7).set(1, I);
        μ.get(7).set(2, 0f);
        μ.get(7).set(3, 0f);
        μ.get(7).set(4, -K);
        μ.get(7).set(5, 0f);
        μ.get(7).set(6, -K);
        μ.get(7).set(7, D);
        μ.get(7).set(8, K);
        μ.get(7).set(9, K);

        μ.get(8).set(0, F);
        μ.get(8).set(1, H);
        μ.get(8).set(2, 0f);
        μ.get(8).set(3, 0f);
        μ.get(8).set(4, -C);
        μ.get(8).set(5, 0f);
        μ.get(8).set(6, -J);
        μ.get(8).set(7, K);
        μ.get(8).set(8, C);
        μ.get(8).set(9, J);

        μ.get(9).set(0, F);
        μ.get(9).set(1, H);
        μ.get(9).set(2, 0f);
        μ.get(9).set(3, 0f);
        μ.get(9).set(4, -J);
        μ.get(9).set(5, 0f);
        μ.get(9).set(6, -C);
        μ.get(9).set(7, K);
        μ.get(9).set(8, J);
        μ.get(9).set(9, C);



    }

    //Matriz tau
    private static void calculatelocalτ(Matrix τ){

        τ.get(0).set(0, 59f);
        τ.get(1).set(0, -1f);
        τ.get(2).set(0, -1f);
        τ.get(3).set(0, -1f);
        τ.get(4).set(0, 4f);
        τ.get(5).set(0, 4f);
        τ.get(6).set(0, 4f);
        τ.get(7).set(0, 4f);
        τ.get(8).set(0, 4f);
        τ.get(9).set(0, 4f);


    }




    //Muestra ls bs
    public static void showbs(ArrayList<Vector> bs){
        for (int i = 0; i < bs.size() ; i++) {
            System.out.println("b del elemento"+(i+1));
            bs.get(i).Show();
            System.out.println("**********************************");
        }
    }

    //Crea el elemento local
    public static float calculateLocalD(int i, Mesh m){
        float D,a,b,c,d,e,f,g,h,i1;

        Element el = m.getElement(i);

        Node n1 = m.getNode(el.getNode1()-1);
        Node n2 = m.getNode(el.getNode2()-1);
        Node n3 = m.getNode(el.getNode3()-1);
        Node n4 = m.getNode(el.getNode4()-1);

        a=n2.getX()-n1.getX();
        b=n2.getY()-n1.getY();
        c=n2.getZ()-n1.getZ();
        d=n3.getX()-n1.getX();
        e=n3.getY()-n1.getY();
        f=n3.getZ()-n1.getZ();
        g=n4.getX()-n1.getX();
        h=n4.getY()-n1.getY();
        i1=n4.getZ()-n1.getZ();
        //SE CALCULA EL DETERMINANTE DE ESTA MATRIZ UTILIZANDO LA REGLA DE SARRUS

        D = a*e*i+d*h*c+g*b*f-g*e*c-a*h*f-d*b*i1;

        return D;
    }

    //Calcula la magnitud de un vector.
    public static float calculateMagnitude(float v1, float v2){
        return (float) Math.sqrt(Math.pow(v1,2)+Math.pow(v2,2));
    }

    public static float calculateLocalArea(int i, Mesh m){
        //Formula de Herón
        float A,s,a,b,c;
        Element e = m.getElement(i);
        Node n1 = m.getNode(e.getNode1()-1);
        Node n2 = m.getNode(e.getNode2()-1);
        Node n3 = m.getNode(e.getNode3()-1);

        a = calculateMagnitude(n2.getX()-n1.getX(),n2.getY()-n1.getY());
        b = calculateMagnitude(n3.getX()-n2.getX(),n3.getY()-n2.getY());
        c = calculateMagnitude(n3.getX()-n1.getX(),n3.getY()-n1.getY());
        s = (a+b+c)/2.0f;

        A = (float) Math.sqrt(s*(s-a)*(s-b)*(s-c));

        return A;
    }
    public static float  ab_ij(float ai, float aj, float al, float bi, float bj, float bl){
        return (ai-al)*(bj-bl)-(aj-al)*(bi-bl);
    }

    //Calcula la la matriz local A
    public static void calculateLocalA(int i,Matrix A, Mesh m){
        Element e = m.getElement(i);
        Node n1 = m.getNode(e.getNode1()-1);
        Node n2 = m.getNode(e.getNode2()-1);
        Node n3 = m.getNode(e.getNode3()-1);
        Node n4 = m.getNode(e.getNode4()-1);

        A.get(0).set(0,ab_ij(n3.getY(),n4.getY(),n1.getY(),n3.getZ(),n4.getZ(),n1.getZ()));
        A.get(0).set(1,ab_ij(n4.getY(),n2.getY(),n1.getY(),n4.getZ(),n2.getZ(),n1.getZ()));
        A.get(0).set(2,ab_ij(n2.getY(),n3.getY(),n1.getY(),n2.getZ(),n3.getZ(),n2.getZ()));
        A.get(1).set(0,ab_ij(n4.getX(),n3.getX(),n1.getX(),n4.getZ(),n3.getZ(),n1.getZ()));
        A.get(1).set(1,ab_ij(n2.getX(),n4.getX(),n1.getX(),n2.getZ(),n4.getZ(),n1.getZ()));
        A.get(1).set(2,ab_ij(n3.getX(),n2.getX(),n1.getX(),n3.getZ(),n2.getZ(),n1.getZ()));
        A.get(2).set(0,ab_ij(n3.getX(),n4.getX(),n1.getX(),n3.getY(),n4.getY(),n1.getY()));
        A.get(2).set(1,ab_ij(n4.getX(),n2.getX(),n1.getX(),n4.getY(),n2.getY(),n1.getY()));
        A.get(2).set(2,ab_ij(n2.getX(),n3.getX(),n1.getX(),n2.getY(),n3.getY(),n1.getY()));




        A.get(0).set(0,n3.getY()-n1.getY());
        A.get(0).set(1, n1.getY()-n2.getY());
        A.get(1).set(0,n1.getX()-n3.getX());
        A.get(1).set(1, n2.getX()-n1.getX());


    }

    //Calcula y llena la matriz B
    private static void calculateB(Matrix B){
        B.get(0).set(0, -1f);
        B.get(0).set(1, 1f);
        B.get(0).set(2, 0f);
        B.get(1).set(3, 0f);
        B.get(1).set(0, -1f);
        B.get(1).set(1, 0f);
        B.get(1).set(2, 1f);
        B.get(1).set(3, 0f);
        B.get(2).set(0, -1f);
        B.get(2).set(1, 0f);
        B.get(2).set(2, 0f);
        B.get(2).set(3, 1f);

    }

   /* public static float calculateLocalVolume(int ind, Mesh m){
        float V,a,b,c,d,e,f,g,h,i;
        Element el = m.getElement(ind);
        Node n1 = m.getNode(el.getNode1()-1);
        Node n2 = m.getNode(el.getNode2()-1);
        Node n3 = m.getNode(el.getNode3()-1);
        Node n4 = m.getNode(el.getNode4()-1);

        a=n2.getX()-n1.getX();
        b=n2.getY()-n1.getY();
        c=n2.getZ()-n1.getZ();
        d=n3.getX()-n1.getX();
        e=n3.getY()-n1.getY();
        f=n3.getZ()-n1.getZ();
        g=n4.getX()-n1.getX();
        h=n4.getY()-n1.getY();
        i=n4.getZ()-n1.getZ();
        //PARA EL DETERMINANTE SE USA LA REGLA  DE SARRUS
        //IMPOSIBLE CONVERTIR DE DE DOUBLE A FLOAT
        //CASTEO A FLOAT
        V= (float) ((1.00/6.00)*(a*e*i+d*h*c+g*b*f-g*e*c-a*h*f-d*b*i));
        return  V;
    }*/

    //Metodo que cera una matriz local K y lo alamcena en m
    private static Matrix createLocalK(int element,Mesh m){
        // K=(k*Ve/D^2)Bd*At*A*B := k_4*4
        Matrix μ = new Matrix();
        calculateLocalμ(CreateC1(element,m),CreateC2(element,m),μ);

        float J;
   //EI es igual a 178

        J = calculateLocalJ(element,m);
        Matrix k = new Matrix(30, 30, 0f);
        for(int i = 0; i<30;i++){
            for(int j = 0; i<30; j++){
                k.get(i).set(j, (178f*J*μ.get(i).get(j)));
                k.get(10+i).set(j+10, (178f*J*μ.get(i).get(j)));
                k.get(20+i).set(j+20, (178f*J*μ.get(i).get(j)));
            }
        }

        return k;
    }

    public static float calculateLocalJ(int i, Mesh m){
        float J,a,b,c,d,e,f,g,h,i1;
        Element el = m.getElement(i);
        Node n1 = m.getNode(el.getNode1()-1);
        Node n2 = m.getNode(el.getNode2()-1);
        Node n3 = m.getNode(el.getNode3()-1);
        Node n4 = m.getNode(el.getNode4()-1);

        a=n2.getX()-n1.getX();
        b=n3.getX()-n1.getX();
        c=n4.getX()-n1.getX();
        d=n2.getY()-n1.getY();
        e=n3.getY()-n1.getY();
        f=n4.getY()-n1.getY();
        g=n2.getZ()-n1.getZ();
        h=n3.getZ()-n1.getZ();
        i1=n4.getZ()-n1.getZ();
        //SE CALCULA EL DETERMIANTE DE ESTA MATRIZ UTILIZANDO METODO DE SARRUS

        J = a*e*i+d*h*c+g*b*f-g*e*c-a*h*f-d*b*i1;

        return J;
    }

    //Metodo que crea el elemento local b y lo almacena en m
    public static Vector createLocalb(int element,Mesh m){
    Matrix τ = new Matrix();
        τ.get(0).set(0, 59f);
        τ.get(1).set(0, -1f);
        τ.get(2).set(0, -1f);
        τ.get(3).set(0, -1f);
        τ.get(4).set(0, 4f);
        τ.get(5).set(0, 4f);
        τ.get(6).set(0, 4f);
        τ.get(7).set(0, 4f);
        τ.get(8).set(0, 4f);
        τ.get(9).set(0, 4f);

        float J;

        J = calculateLocalJ(element,m);
        Matrix ma = new Matrix(30, 1, 0f);
        //f es corrdenada en x=-15  y=77  z=28
        Vector f = new Vector();

        f.add(0,-15f);
        f.add(0,77f);
        f.add(0,28f);


        for (int i =0; i<30;i++){
            ma.get(i).set(0, (J/120f*τ.get(i).get(0)));
            ma.get(10+i).set(1, (J/120f*τ.get(i).get(0)));
            ma.get(20+i).set(2, (J/120f*τ.get(i).get(0)));
        }

        Vector b = new Vector(30,0f);
        productMatrixVector(ma,f,b);


        return b;
    }

    //Esta funcion crea los sitemas locales (K y b) y almacena los datos en sus respectivas listas
    public static void crearSistemasLocales(Mesh m, ArrayList<Matrix> localKs, ArrayList<Vector> localbs){

        for(int i = 0; i<m.getSize(Sizes.ELEMENTS.ordinal()); i++){
            localKs.add(createLocalK(i,m));
            localbs.add(createLocalb(i,m));
        }
    }

    //Esta funcion realiza el ensamblaje de la matriz K global, recibe el elemento actual, la matriz K local
    //y la matriz K global en la cual se realizara el ensamblaje
    public static void assemblyK(Element e, Matrix localK, Matrix K){
        int index1 = e.getNode1() - 1;
        int index2 = e.getNode2() - 1;
        int index3 = e.getNode3() - 1;
        int index4 = e.getNode4() - 1;

        K.get(index1).set(index1, K.get(index1).get(index1) + localK.get(0).get(0));
        K.get(index1).set(index2, K.get(index1).get(index2) + localK.get(0).get(1));
        K.get(index1).set(index3, K.get(index1).get(index3) + localK.get(0).get(2));
        K.get(index1).set(index4, K.get(index1).get(index4) + localK.get(0).get(3));
        K.get(index2).set(index1, K.get(index2).get(index1) + localK.get(1).get(0));
        K.get(index2).set(index2, K.get(index2).get(index2) + localK.get(1).get(1));
        K.get(index2).set(index3, K.get(index2).get(index3) + localK.get(1).get(2));
        K.get(index2).set(index4, K.get(index2).get(index4) + localK.get(1).get(3));
        K.get(index3).set(index1, K.get(index3).get(index1) + localK.get(2).get(0));
        K.get(index3).set(index2, K.get(index3).get(index2) + localK.get(2).get(1));
        K.get(index3).set(index3, K.get(index3).get(index3) + localK.get(2).get(2));
        K.get(index3).set(index4, K.get(index3).get(index4) + localK.get(2).get(3));
        K.get(index4).set(index1, K.get(index4).get(index1) + localK.get(3).get(0));
        K.get(index4).set(index2, K.get(index4).get(index2) + localK.get(3).get(1));
        K.get(index4).set(index3, K.get(index4).get(index3) + localK.get(3).get(2));
        K.get(index4).set(index4, K.get(index4).get(index4) + localK.get(2).get(3));

    }

    //Esta funcion realiza el ensamblaje del arreglo b global, recibe el elemento actual, el arreglo b local
    //y el arreglo b glocal en el cual se realizara el ensablaje
    public static void assemblyb(Element e, Vector localb, Vector b){
        int index1 = e.getNode1() - 1;
        int index2 = e.getNode2() - 1;
        int index3 = e.getNode3() - 1;
        int index4 = e.getNode4() - 1;

        b.set(index1, b.get(index1) + localb.get(0));
        b.set(index2, b.get(index2) + localb.get(1));
        b.set(index3, b.get(index3) + localb.get(2));
        b.set(index4, b.get(index4) + localb.get(3));
    }

    //Se realiza el ensamblaje de los sistemas locales K y B utilizando las funciones assemblyK y assemblyb
    public static void ensamblaje(Mesh m, ArrayList<Matrix> localKs, ArrayList<Vector> localbs, Matrix K,Vector b){
        for(int i=0; i<m.getSize(Sizes.ELEMENTS.ordinal()); i++){
            Element e = m.getElement(i);
            assemblyK(e,localKs.get(i),K);
            assemblyb(e,localbs.get(i),b);
        }
    }

    //Funcion que aplica la condicion de neumann al vector b
    public static void applyNeumann(Mesh m,Vector b){
        for(int i=0;i <m.getSize(Sizes.NEUMANN.ordinal()); i++){
            Condition c = m.getCondition(i,Sizes.NEUMANN);
            b.set(c.getNode1()-1, b.get(c.getNode1()-1) + c.getValue());
        }
    }

    //Funcion que aplica la condicion de dirichlet al sistema de ecuaciones
    public static void applyDirichlet(Mesh m,Matrix K,Vector b){
        for(int i=0; i<m.getSize(Sizes.DIRICHLET.ordinal()); i++){

            Condition c = m.getCondition(i,Sizes.DIRICHLET);
            int index = c.getNode1()-1;

            K.remove(index);
            b.remove(index);

            for(int row=0; row < K.size(); row++){
                float cell = K.get(row).get(index);
                K.get(row).remove(index);
                b.set(row, b.get(row) + (-1*c.getValue()) * cell);
            }
        }
    }

    //Funcion que calcula el resultado del SEL
    public static void calculate(Matrix K, Vector b, Vector T){
        System.out.println("Iniciando calculo de respuesta...");
        Matrix Kinv = new Matrix();
        System.out.println("Calculo de inversa...");
        inverseMatrix(K, Kinv);
        System.out.println("Calculo de respuesta...");
        productMatrixVector(Kinv, b, T);
    }

}
