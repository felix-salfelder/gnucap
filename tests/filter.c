#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#define untested() (fprintf(stderr, "@@# untested \n@@@:%s:%u:%s\n", \
			   __FILE__, __LINE__, __func__))

int main(int argc, char * const * argv)
{
	if(argc<2){
		exit(1);
	}
	int buf[10];

	FILE * datafile = fopen(argv[1],"r");

	int t = getchar();
	int data = fgetc(datafile);
	int s = 0;
	int e = 0;
	int end = ' ';
	while( (data != EOF) && (t != EOF) ){
		if(e) {
			if (e>9){
				break;
			} else if (t!=' '){
				buf[e] = t;
				e++;
			}
			if (s) { // compare buf
				if (s>=e) {
					e = s = 0;
				} else if (data==' ') {
					e = s = 0;
				} else if (s<10 && data==buf[s]){
					data = t;
					s++;
				} else {

				}
			} else { // copy to data
				if (data==buf[0]) {
					s = 1;
				}
				data = t;
			}
		}else if (t=='~' && '0' <= data && data <= '9'){
			data = '~';
			s = 1;
		} else if(t=='~' && ( 'u' == data || 'n' == data || 'f' == data || 'p' == data || 'K' == data )){
			end = data;
			data = '~';
			s = 1;
		} else if (s) {
			if (t==' ' && '0' <= data && data <= '9') {
				data = ' ';
				end=' ';
			} else if (t=='E' && '0' <= data && data <= '9' ) {
				data = 'E';
				buf[0] = 'E';
				e = 1;
				s = 0;
			} else if (t==end){
				data = end;
				end = ' ';
				s = 0;
			} else if ((t=='n' || t=='f' || t=='p' || t=='u' || t=='K') && '0' <= data && data <= '9') {
				data = t;
				end = ' ';
				s = 1;
			}
		}else if(t=='~' && (data == '-' || data == ' ')){
			data = '~';
		}else{
			assert(!s);
		}
		putchar(data);
	   data = fgetc(datafile);
		t = getchar();
		if(data == '\n'){
			while(t!='\n' && t!=EOF){
				t = getchar();
			}
		}else if(t == '\n'){
			while(data!='\n' && data!=EOF){
				putchar(data);
				data = fgetc(datafile);
			}
		}
		if(data=='\n' || t=='\n'){
			s = 0;
			e = 0;
			end = ' ';
		}
	}
	while( data != EOF ){
		putchar(data);
	   data = fgetc(datafile);
	}
	fclose(datafile);
	return 0;
}
