#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define SIZE 20

typedef struct{
	float cofactor;
	int degree;
}term;
void bisection();
void regulaFalsi();
void newtonRaphson();
void inverseMatrix();	
void gaussJordan();
void gaussSeidel();
void numericDifferantial();
void simpson();
void trapez();
void gregoryNewton();
void choice();
void addRow(float matrix[SIZE][SIZE],int row, int addingRow, int length, float addingMultiplier);
void multiplyRow(float matrix[SIZE][SIZE],int row, int length, float multiplier);
void takePolynomial(term polynomial[], int *length);
void takeRange(float *start, float *end);
void differantial(term polynomial[], term diffpolynomial[], int length);
void adjoint(float matrix[SIZE][SIZE], int length, float adjMatrix[SIZE][SIZE]);
void printMatrix(float matrix[SIZE][SIZE], int length);
void inverse(float adjMatrix[SIZE][SIZE], float invMatrix[SIZE][SIZE], float det, int length);
float absDiff(float start, float end);
float functionValue(float input, int length, term polynomial[]);
float determinant(float matrix[SIZE][SIZE], int length);
int factorial(int a);

int main(int argc, char *argv[]) {
	choice();
	return 0;
}

void choice(){
	int n;
	printf("\n\t**********NUMERICAL ANALYSIS PROJECT**********\n\n\n\t1. Bisection method\n\t2. Regula-Falsi method\n\t3. Newton-Raphson method\n\t4. Inverse of NxN matrix\n\t5. Gauss Elemination method\n\t6. Gauss Seidal Method\n\t7. Numerical Differantiation (higher, lower and symmetric difference)\n\t8. Simpson method\n\t9. Trapez method\n\t10. Gregory Newton Enterpolation\n\n\tType the number of method that you want to choose, Type 0 for quit : ");	
	scanf("%d",&n);
	switch (n){
		case 1:
			bisection();
			break;
		case 2:
			regulaFalsi();
			break;
		case 3:
			newtonRaphson();
			break;
		case 4:
			inverseMatrix();
			break;
		case 5:
			gaussJordan();
			break;
		case 6:
			gaussSeidel();
			break;
		case 7:
			numericDifferantial();
			break;
		case 8: 
			simpson();
			break;
		case 9:
			trapez();
			break;
		case 10: 
			gregoryNewton();
			break;
		case 0:
			exit(1);
	}
}

void takePolynomial(term polynomial[], int *length){
	int i;
	system("cls");
	printf("\n\tEnter term number of polynomial : ");
	scanf("%d",&*length);
	for(i=0; i<*length; i++){
		printf("\tEnter the cofactor of the %d. term : ",i+1);
		scanf("%f",&polynomial[i].cofactor);
		printf("\tEnter the power of the %d. term : ",i+1);
		scanf("%d",&polynomial[i].degree);
	}
	system("cls");
	printf("\n\tThis is the polynomial you entered : ");
	for(i=0; i<*length; i++){
		printf("%.2fx^%d + ",polynomial[i].cofactor,polynomial[i].degree);
	}
}

float functionValue(float input, int length, term polynomial[]){
	int i; 
	float sum = 0;
	for(i=0; i<length; i++){
		sum+= polynomial[i].cofactor*pow(input,polynomial[i].degree); 
	}
	return sum;
}

void bisection(){
	int i=1, length;
	float start, end, middle, epsilon;
	term polynomial[SIZE];
	
	takePolynomial(polynomial,&length);
	takeRange(&start,&end);
	while(functionValue(start,length,polynomial)*functionValue(end,length,polynomial)>0){
		printf("\n\n\tThere is no real root in the value range you entered, please try a different range.");
		takeRange(&start,&end);
	}
	printf("\n\tThis is the polynomial you entered : ");
	for(i=0; i<length; i++){
		printf("%.2fx^%d + ",polynomial[i].cofactor,polynomial[i].degree);
	}
	printf("\n\n\tEnter the epsilon(margin of error) : ");
	scanf("%f",&epsilon);
	while(absDiff(start,end)/pow(2,i)>epsilon){
		middle = (start+end)/2;
		if(functionValue(start,length,polynomial)*functionValue(middle,length,polynomial)<0){
			end = middle;
		}
		else{
			start = middle;
		}
		printf("\n\tRange processed in the %d. iteration : [(%.2f,%.2f) , (%.2f,%.2f)]",i,start,functionValue(start,length,polynomial),end,functionValue(end,length,polynomial));
		i++;	
	}
	printf("\n\n\tThe value obtained as a approximately result of the %d. iteration = %.2f",i-1,(start+end)/2);		
}

void takeRange(float *start, float *end){
	printf("\n\n\tEnter the start value of the range you want to process : ");
	scanf("%f",&*start);
	printf("\tEnter the end value of the range you want to process : ");
	scanf("%f",&*end);
	system("cls");
	printf("\n\tThis is the range you entered : (%.2f,%.2f)",*start,*end);
}

float absDiff(float start, float end){
	if(start-end<0){
		return end-start;
	}
	else{
		return start-end;
	}
}

void regulaFalsi(){
	int i=1, length;
	float start, end, middle, epsilon;
	term polynomial[SIZE];
	
	takePolynomial(polynomial,&length);
	takeRange(&start,&end);
	while(functionValue(start,length,polynomial)*functionValue(end,length,polynomial)>0){
		printf("\n\n\tThere is no real root in the value range you entered, please try a different range.");
		takeRange(&start,&end);
	}
	printf("\n\n\tEnter the epsilon(margin of error) : ");
	scanf("%f",&epsilon);
	printf("\n\tThis is the polynomial you entered : ");
	for(i=0; i<length; i++){
		printf("%.2fx^%d + ",polynomial[i].cofactor,polynomial[i].degree);
	}
	while(absDiff(start,end)/pow(2,i)>epsilon){
		middle = (end*functionValue(start,length,polynomial)-(start*functionValue(end,length,polynomial)))/(functionValue(start,length,polynomial)-functionValue(end,length,polynomial));
		if(functionValue(start,length,polynomial)*functionValue(middle,length,polynomial)<0){
			end = middle;
		}
		else{
			start = middle;
		}
		printf("\n\tRange processed in the %d. iteration : [(%.2f,%.2f) , (%.2f,%.2f)]",i,start,functionValue(start,length,polynomial),end,functionValue(end,length,polynomial));
		i++;	
	}
	printf("\n\n\tThe value obtained as a approximately result of the %d. iteration = %.2f",i-1,(end*functionValue(start,length,polynomial)-(start*functionValue(end,length,polynomial)))/(functionValue(start,length,polynomial)-functionValue(end,length,polynomial)));		
}

void differantial(term polynomial[], term diffpolynomial[], int length){
	int i;
	for (i=0; i<length; i++){
		diffpolynomial[i].degree = polynomial[i].degree-1;
		diffpolynomial[i].cofactor = polynomial[i].cofactor*polynomial[i].degree;
	}
}

void newtonRaphson(){
	int i=1, length;
	float start, next, epsilon;
	term polynomial[SIZE], diffpolynomial[SIZE];
	
	takePolynomial(polynomial,&length);
	differantial(polynomial,diffpolynomial,length);
	printf("\n\n\tBaşlangıç değerini giriniz : ");
	scanf("%f",&start);
	printf("\n\tEnter the epsilon(margin of error) : ");
	scanf("%f",&epsilon);
	printf("\n\tThis is the polynomial you entered : ");
	for(i=0; i<length; i++){
		printf("%.2fx^%d + ",polynomial[i].cofactor,polynomial[i].degree);
	}
	next = start - (functionValue(start,length,polynomial)/functionValue(start,length,diffpolynomial));
	while(absDiff(start,next)>epsilon){
		printf("\n\tThe value obtained as a result of the %d. iteration : %f",i,next);			
		start = next;
		next = next - (functionValue(next,length,polynomial)/functionValue(next,length,diffpolynomial));
		i++;
		if(i==100){
			printf("\n\tA valid value could not be found in 100 iterations with the Regula Falsi method.");
			break;
		}	
	}
	if(i<100){
	printf("\n\tThe value obtained as a result of the %d. iteration : %f",i,next);
	}
}

void numericDifferantial(){
	int length, choice;
	float epsilon, point, result;
	term polynomial[SIZE];
	
	takePolynomial(polynomial,&length);
	printf("\n\tEnter the x value of that you want to find  : ");
	scanf("%f",&point);
	printf("\tEnter the epsilon(margin of error) : ");
	scanf("%f",&epsilon);
	printf("\n\tType 1 for higher difference,\n\tType 2 for lower difference,\n\tType 3 for symmetric difference : ");
	scanf("%d",&choice);
	if(choice == 1){
		result = (functionValue(point+epsilon,length,polynomial)-functionValue(point,length,polynomial))/epsilon;
	}
	else if(choice == 2){
		result = (functionValue(point,length,polynomial)-functionValue(point-epsilon,length,polynomial))/epsilon;
	}
	else if(choice == 3){
		result = (functionValue(point+length,length,polynomial)-functionValue(point-epsilon,length,polynomial))/epsilon*2;
	}
	printf("Differantial value of \t%f point : %f",point,result);
}

void trapez(){
	int i, length, n;
	float start, end, middle, epsilon, h, result, sum = 0;
	term polynomial[SIZE];
	
	takePolynomial(polynomial,&length);
	takeRange(&start,&end);
	printf("\n\tEnter number of divisions : ");
	scanf("%d",&n);
	h = (end-start)/n;
	for(i=1; i<n; i++){
		sum += functionValue(start+(h*i),length,polynomial);
	}
	result = h*((functionValue(start,length,polynomial)+functionValue(start,length,polynomial))/2+sum);
	printf("\n\tValue of integral : %f",result);
}

void simpson(){
	int i, length, n, cho;
	float start, end, middle, epsilon, h, result, sum = 0;
	term polynomial[SIZE];
	
	takePolynomial(polynomial,&length);
	takeRange(&start,&end);
	printf("\n\tType 1 for 1/3 rule, type 2 for 3/8 rule : ");
	scanf("%d",&cho);
	if(cho == 1){
		printf("\n\tEnter number of divisions : ");
		scanf("%d",&n);
		h = (end-start)/n;
		for(i=1; i<n; i++){
			if(i%2==1){
				sum += 4*(functionValue(start+(h*i),length,polynomial));
			}
			else{
				sum += 2*(functionValue(start+(h*i),length,polynomial));
			}
		}
		result = (h/3)*(functionValue(start,length,polynomial)+functionValue(end,length,polynomial)+sum);
	}
	else if(cho == 2){
		h = (end-start)/3;
		result = (end-start) * (functionValue(start,length,polynomial) + 3*functionValue(start+h,length,polynomial) + 3*functionValue(start+2*h,length,polynomial) + functionValue(end,length,polynomial)) / 8;
	}
	printf("\n\tValue of integral : %f",result);
}

void inverseMatrix(){
	float matrix[SIZE][SIZE], adjMatrix[SIZE][SIZE], invMatrix[SIZE][SIZE];
	int length, i, j, sum = 0;
	float det;
	
	system("cls");
	printf("\n\tEnter the N value of NxN matrix : ");
	scanf("%d",&length);
		
	for(i=0; i<length; i++){
		for(j=0; j<length; j++){
			printf("\tEnter the element [%d,%d] : ",i+1,j+1);
			scanf("%f",&matrix[i][j]);
		}
	}
	
	system("cls");
	printf("\n\tThis is your matrix :\n");
	adjoint(matrix,length,adjMatrix);
	printMatrix(matrix,length);
	if(determinant(matrix, length) == 0){
		printf("\n\tInverse of a matrix does not exists if its determinant is equal to 0.");
	}
	else{
		printf("\n\tMatrisin determinant değeri : %f\n\n\tKatsayılar matrisi :\n",determinant(matrix, length));
		printMatrix(adjMatrix,length);
		printf("\n\tGirilen matrisin tersi :\n");
		inverse(adjMatrix, invMatrix, determinant(matrix, length), length);
		printMatrix(invMatrix, length);
	}
	}	

void gaussJordan(){
	float matrix[SIZE][SIZE], values[SIZE];
	float sum;
	int i, j, length;
	
	system("cls");
	printf("\n\tEnter the N value of NxN matrix : ");
	scanf("%d",&length);
	for(i=0; i<length; i++){
		for(j=0; j<length+1; j++){
			printf("\tEnter the element [%d,%d] : ",i+1,j+1);
			scanf("%f",&matrix[i][j]);
		}
	}
	for(i=0; i<length; i++){
		multiplyRow(matrix,i,length,1/matrix[i][i]);
		for(j=i+1; j<length; j++){
			addRow(matrix,j,i,length,-(matrix[j][i]/matrix[i][i]));
		}
	}
	values[length-1] = matrix[length-1][length];
	for(i=length-2; i>=0; i--){
		sum = 0;
		for(j=length-1; j>i; j--){
			sum -= matrix[i][j]*values[j];
		}
		sum += matrix[i][length];
		values[i] = sum;
	}
	
	system("cls");
	printf("\n\tThis echelon form of your matrix :");
	printMatrix(matrix,length);
	printf("\n\tThis is solvation of your matrix :");
	for(i=0; i<length; i++){
		printf("\tx%d = %f",i+1,values[i]);
	}
}

void gaussSeidel(){
	float sum = 0, epsilon, error = 0.1;
	float matrix[SIZE][SIZE];
	float values[SIZE];
	int i, j, k, length, count = 1;
	
	system("cls");
	printf("\n\tEnter the N value of NxN matrix : ");
	scanf("%d",&length);
		
	for(i=0; i<length; i++){
		for(j=0; j<length+1; j++){
			printf("\tEnter the element [%d,%d] : ",i+1,j+1);
			scanf("%f",&matrix[i][j]);
		}
	}
	
	for(i=0; i<length; i++){
		printf("Enter the first value of \t%d. variable : ",i+1);
		scanf("%f",&values[i]);
	}
	
	printf("\n\n\tEnter the epsilon(margin of error) : ");
	scanf("%f",&epsilon);
	system("cls");
	while((error>epsilon)&&(count<50)){
		for(i=0; i<length; i++){
			for(j=0; j<length; j++){
				if(i!=j){
					sum -= matrix[i][j]*values[j];
				}
			}
			sum += matrix[i][length];
			error = absDiff(sum/matrix[i][i],values[i]);
			if(error<absDiff(sum/matrix[i][i],values[i])){
				error = absDiff(sum/matrix[i][i],values[i]);
			}
			values[i] = sum/matrix[i][i];
			sum = 0;
		}
		printf("\n\t The values obtained as a result of the %d. iteration :",count);
		for(i=0; i<length; i++){
		printf("\tx%d = %f",i+1,values[i]);
		}
		count++;
	}
	if(count == 50){
		printf("\n\tValues diverged.");
	}
}

void gregoryNewton(){
	float x[SIZE], y[SIZE][SIZE];
	float h, desiredX;
	int length, i, j;
	term polynomial[SIZE];

	system("cls");
	printf("\n\tEnter h that is the difference of x values : ");
	scanf("%f",&h);
	printf("\tEnter number of values : ");
	scanf("%d",&length);
	printf("\tEnter the x value that you want to find : ");
	scanf("%f",&desiredX);
	printf("\tEnter x value of starting point : ");
	scanf("%f",&x[0]);
	printf("\tEnter y value of starting point : ");
	scanf("%f",&y[0][0]);
	for(i=1; i<length; i++){
		x[i] = x[0]+i*h;
		printf("\tEnter y value of %f point : ",x[i]);
		scanf("%f",&y[0][i]);
	}
	system("cls");
	for(i=0; i<4; i++){
		for(j=0; j<length-i-1; j++){
			y[i+1][j] = y[i][j+1] - y[i][j];
		}
	}
	for(i=0; i<5; i++){
		printf("\nDelta^%df(x)",i);
		for(j=0; j<length-i; j++){
			printf("\t%f",y[i][j]);
		}
	}
	polynomial[0].degree = 4;
	polynomial[0].cofactor = y[4][0]/(pow(h,4)*factorial(4));
	polynomial[1].degree = 3;
	polynomial[1].cofactor = (y[4][0]/(pow(h,4)*factorial(4)))*(-x[0]-x[1]-x[2]-x[3])+(y[3][0]/(pow(h,3)*factorial(3)));
	polynomial[2].degree = 2;
	polynomial[2].cofactor = (y[4][0]/(pow(h,4)*factorial(4)))*(-x[3]*(x[0]+x[1]+x[2])+x[0]+x[1]+x[0]*x[1])+(y[3][0]/(pow(h,3)*factorial(3))*(-x[0]-x[1]-x[2]))+(y[2][0]/(pow(h,2)*factorial(2)));
	polynomial[3].degree = 1;
	polynomial[3].cofactor = (y[4][0]/(pow(h,4)*factorial(4)))*(x[3]*(x[0]+x[1]+x[0]*x[1])-x[0]*x[1]*x[2])+(y[3][0]/(pow(h,3)*factorial(3))*(x[0]+x[1]+x[0]*x[1]))+(y[2][0]/(pow(h,2)*factorial(2)))*(-x[0]-x[1])+(y[1][0]/pow(h,1));
	polynomial[4].degree = 0;
	polynomial[4].cofactor = (y[4][0]/(pow(h,4)*factorial(4)))*(x[0]*x[1]*x[2]*x[3])+(y[3][0]/(pow(h,3)*factorial(3))*(-x[0]*x[1]*x[2])+(y[2][0]/(pow(h,2)*factorial(2)))*(x[0]*x[1]))+(y[1][0]/pow(h,1))*(-x[0])+(y[0][0]);
	printf("\n\tpolynomial : ");
	for(i=0; i<5; i++){
		printf("%.2fx^%d + ",polynomial[i].cofactor,polynomial[i].degree);
	}
	printf("\n\tf(%f) value : %f",desiredX,functionValue(desiredX,5,polynomial));
}

float determinant(float matrix[SIZE][SIZE], int length){
	int i, j, k, a, b;
	float sum = 0;
	
	if(length == 2){
		return (matrix[0][0]*matrix[1][1])-(matrix[0][1]*matrix[1][0]);
	}
	else{
		for(i=0; i<length; i++){
			float minorMatrix[SIZE][SIZE];
			a = 0;
			for(j=1; j<length; j++){
				for(k=0; k<length; k++){
					if(k!=i){
						minorMatrix[a][b] = matrix[j][k];
						b++;
					}
				}
				a++;
				b = 0;
			}
			sum += pow(-1,i)*(matrix[0][i])*(determinant(minorMatrix, length-1));
		}
		return sum;
	}
}

void adjoint(float matrix[SIZE][SIZE], int length, float adjMatrix[SIZE][SIZE]){
	float tmpMatrix[SIZE][SIZE];
	int i, j, k, t, a, b;
	float tmp;
	for(i=0; i<length; i++){
		for(j=0; j<length; j++){
			a = 0;
			b = 0;
			for(k=0; k<length; k++){
				b=0;
				for(t=0; t<length; t++){
					if((k!=i)&&(t!=j)){
						tmpMatrix[a][b] = matrix[k][t];
						if (b < (length - 2)){
				             b++;
				    	}
				        else{
				               b = 0;
				               a++;
				        }
					}
				}				
			}
			adjMatrix[j][i] = determinant(tmpMatrix,length-1);
		}
		
	}
}

void printMatrix(float matrix[SIZE][SIZE], int length){
	int i,j;
	printf("\n");
	for(i=0; i<length; i++){
		for(j=0; j<length; j++){
			printf("\t%f",matrix[i][j]);
			}
		printf("\n");
		}
}

void inverse(float adjMatrix[SIZE][SIZE], float invMatrix[SIZE][SIZE], float det, int length){
	int i,j;
	for(i=0; i<length; i++){
		for(j=0; j<length; j++){
			invMatrix[i][j] = adjMatrix[i][j]/det;
			if((i+j)%2==1){
				invMatrix[i][j] *= -1;
			}
		}
	}
}

void multiplyRow(float matrix[SIZE][SIZE],int row, int length, float multiplier){
	int i;
	for(i=0; i<length+1; i++){
		matrix[row][i] *= multiplier;
	}
}

void addRow(float matrix[SIZE][SIZE],int row, int addingRow, int length, float addingMultiplier){
	int i;
	for(i=0; i<length+1; i++){
		matrix[row][i] += addingMultiplier*matrix[addingRow][i];
	}
}

int factorial(int a){
	if(a<2){
		return 1;
	}
	else{
		return a*factorial(a-1);
	}
}
