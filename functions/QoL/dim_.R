dim_<-function(input){
  #better version of dim, less R nonsense
  if(is.list(input)){
    print(paste0("This list contains: ", str(metaData)))
  }else if(is.matrix(input)){
    print(paste0("This matrix has dimensions: ", dim(input)[1],' rows, and ' ,dim(input)[2], ' columns.' ))
  }else if(is.vector(input)){
    print(paste0("This vector has length: " , length(input)))
  }else{
    print(paste0("why did you check the dimension of a ", typeof(input), " object?"))
  }
}
