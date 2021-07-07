package Classes;

public class Node extends Item{


    //Metodo que nos ayudara a crear las listas de nodos
    public static Node[] createNodes(int n){
        Node[] list = new Node[n];
        for (int i = 0; i < n; i++) {
            list[i] = new Node();
        }
        return list;
    }
    //la cantidad de parametros que se le manda a set value
    @Override
    public void  setValues(int a,float b,float c,float d,int e,int f,int g, int h, float i) {
        id = a;
        x = b;
        y = c;
        //agtregado
        z = d;
    }
}
