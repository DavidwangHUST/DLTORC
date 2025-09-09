template<typename T>
int parseNumber(FILE* input, const char* keyword, const size_t bufLen, T* n){

	char buf[bufLen];
	char buf1[bufLen];
	char comment[1];
	char *ret;

	while (fgets(buf,bufLen, input)!=NULL){
		comment[0]=buf[0];
		if(strncmp(comment,"#",1)==0){
			//printf("Comment!:%s\n",buf);
		}
		else{
			ret=strtok(buf,"=");
			if(strcmp(ret,keyword)==0){
				/*offset buf by keyword size + 1 for the "="*/
				strncpy (buf1, buf+strlen(keyword)+1, bufLen);	
       				printf("%30s: ",keyword);
				getFromString(buf1,n);
				rewind(input);
				return(0);
			}
		}
	}
	rewind(input);
	return(-1);
}

template<typename T>
int parseArray(FILE* input, const char* keyword, const size_t bufLen,
		const size_t arrLen, T y[]){

	char buf[bufLen];
	char buf1[bufLen];
	char comment[1];
	char *ret;

	while (fgets(buf,bufLen, input)!=NULL){
		comment[0]=buf[0];
		if(strncmp(comment,"#",1)==0){
			//printf("Comment!:%s\n",buf);
		}
		else{
			ret=strtok(buf,"=");

			if(strcmp(ret,keyword)==0){
				/*offset buf by keyword size + 1 for the "="*/
				strncpy (buf1, buf+strlen(keyword)+1, bufLen);	
       				printf("%30s:\n",keyword);
				ret=strtok(buf1,",");
				size_t j=0;
				while(ret!=NULL){
					if(j<arrLen){
						//y[j]=atof(ret);
						getFromString(ret,&y[j]);
					}
   		  			ret=strtok(NULL,",");
					j++;
				}
				rewind(input);
				if(j!=arrLen){
					printf("Check no: of values entered for %s\n",keyword);
					printf("%lu values required!\n",arrLen);
					return(-1);
				}
				else{
					return(0);
				}
			}
		}
	}
	rewind(input);
	return(-1);
}

