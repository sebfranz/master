dim_<-function(input){
  #better version of dim, less R nonsense
  if(is.list(input)){
    str(metaData)
  }else if(is.matrix(input)){
    dim(input)
  }else if(is.vector(input)){
    length(input)
  }else{
    print(paste0("why did you check the dimension of a ", typeof(input), "object?"))
  }
}
