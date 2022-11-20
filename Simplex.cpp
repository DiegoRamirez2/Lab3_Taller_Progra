#include "Simplex.h"


// estas funciones deben sacarlas de aqui
// todas las funciones externas deben salir
//extern "C" int *ivector(long nl, long nh);
extern "C" float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch);
extern "C" void simplx(float **a, int m, int n, int m1, int m2, int m3, int *icase,int izrov[], int iposv[]);

Simplex::Simplex(int m1, int m2, int m3){
	// faltan cosas aqui
	this->solved=false;
}

Simplex::~Simplex(){

}

bool Simplex::loadTxt(string filename){
	string line;
	// neceitamos la libreria de manejo de archivo fstream
	// sstream stream de string
	fstream txt_stream(filename);
	stringstream charStream;
	string number_txt;
	vector<float> *number_float;
	
	if (!txt_stream.is_open())	{
		cout<<"Archivo: " << filename << " no encontrado" << endl;
		return false;
	} else {
		cout<<"Archivo encontrado " << endl;
	}
	// notar que getline retorna verdadero cuando logra leer algo
	getline(txt_stream, line); // lee caracteres hasta fin de linea
	//cout<< line << endl;
	charStream <<line;
	getline(charStream, number_txt, ' '); // idem hasta espacio
	m1=stoi(number_txt); // transforma string a int
	getline(charStream, number_txt, ' '); // idem hasta espacio
	m2=stoi(number_txt); // transforma string a int
	getline(charStream, number_txt, ' '); // idem hasta espacio
	m3=stoi(number_txt); // transforma string a int

	// lee mientras exista una linea que leer
	while(getline(txt_stream,line)) {
		charStream.clear();
		charStream << line;
		number_float=new vector<float>;
		while(getline(charStream, number_txt, ' ')) {
			number_float->push_back(stof(number_txt));
		}
		a.push_back(*number_float);
	}
	N=a[0].size()-1;
	M=a.size()-1;
	NP= N+1;
	MP=M+2;
	NM1M2 = N+m1+m2;

	// falta chequear que sea arreglo rectangular 
	// ie todas las filas del mismo tamaño
	// m1+m2+m3 == M

	// ultimo detalle agregar fila con ceros
	number_float=new vector<float>;
	for(int i=0;i<=N; i++){
		number_float->push_back(0.0);
	}
	a.push_back(*number_float);

	txt_stream.close();
    return true;
}

void Simplex::print(){
	cout<<"m1: "<<m1<< " m2:" << m2 << " m3:" << m3<< endl;
	cout<<"N:"<<N<< " M:"<<M<<endl;
	cout<<"NP:"<<NP<< " MP:"<<MP<<endl;
	
	int i=0, j=0;
	for(i=0; i<a.size(); i++){
		for(j=0;j<a[i].size();j++){
			cout<<a[i][j] <<" ";
		}
		cout<<endl;
	}
}

void Simplex::insertConstraint(int b, int var, int type){

}

Simplex* Simplex::clone() {
    Simplex *s1 = new Simplex(m1, m2, m3);
	// hacer la copia
	return s1;
}

void Simplex::solve() {
/*
#define NP 5        // NP >= N+1 
#define MP 6        // MP >= M+2 
*/
	int i,icase,j,izrov[N+1],iposv[M+1]; // esto lo cambie izrov e iposv
	// por el momento no cargamos nada del txt y resolvemos este sistema
	float c[MP][NP]=
		{0.0,1.0,1.0,3.0,-0.5,
		740.0,-1.0,0.0,-2.0,0.0,
		0.0,0.0,-2.0,0.0,7.0,
		0.5,0.0,-1.0,1.0,-2.0,
		9.0,-1.0,-1.0,-1.0,-1.0,
		0.0,0.0,0.0,0.0,0.0};
	float **a;

	//static char *txt[NM1M2+1]=
	// OJO: generar un arreglo de string con NM1M2+1 elementos 
	string txt[NM1M2+1]=
				{" ","x1","x2","x3","x4","y1","y2","y3"};

	//izrov=ivector(1,N);
	//iposv=ivector(1,M);
	a=convert_matrix(&c[0][0],1,MP,1,NP);
/*    a=(float **) calloc(MP+1,sizeof(float*)); 
    for(int i=0;i<=M;i++){
        a[i]=(float *) calloc(NP+1,sizeof(float));
        for(j=1;j<=MP;++j){
            a[i][j]=c[i+1][j-1];
        }
    }
    */
	simplx(a,M,N,m1,m2,m3,&icase,izrov,iposv);

	// dejar todo esto en la funcion print y cambiar printf por cout << ... 
	if (icase == 1)
		printf("\nunbounded objective function\n");
	else if (icase == -1)
		printf("\nno solutions satisfy constraints given\n");
	else {
		printf("\n%11s"," ");
		for (i=1;i<=N;i++)
			if (izrov[i] <= NM1M2) 
				cout << txt[izrov[i]];
		printf("\n");
		for (i=1;i<=M+1;i++) {
			if (i == 1 || iposv[i-1] <= NM1M2) {
				if (i > 1)
					cout << txt[iposv[i-1]];
				else
					printf("  ");
				printf("%10.10f",a[i][1]);
				for (j=2;j<=N+1;j++)
					if (izrov[j-1] <= NM1M2)
						printf("%10.2f",a[i][j]);
				printf("\n");
			}
		}
	}
	// deben hacer el destructor

	//free_convert_matrix(a,1,MP,1,NP);
	//free_ivector(iposv,1,M);
	//free_ivector(izrov,1,N);
	
}

void Simplex::simplx_tmp() {

}

void Simplex::simp1_tmp() {

}
void Simplex::simp2_tmp() {

}
void Simplex::simp3_tmp() {

}
