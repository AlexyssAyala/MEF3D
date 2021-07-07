package Classes;

public class Condition extends Item{




    //Metodo que nos ayudara a crear las listas de condiciones
    public static Condition[] createConditions(int n){
        Condition[] list = new Condition[n];
        for (int i = 0; i < n; i++) {
            list[i] = new Condition();
        }
        return list;
    }

    //la cantidad de parametros que se le manda a set value
    @Override
    public void setValues(int a,float b,float c,float d,int e,int f,int g, int h, float i) {
        node1 = e;
        value = i;
    }



}
