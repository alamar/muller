
// macro definition for referring to array in class like o.a[i] from python
// define %array_class(smth, smthArray); first !
// usage example: ARRAYMEMBER(Spam, array, floatArray);

%define ARRAYMEMBER(cls, name, type)

%extend cls{
    type * name ## _getarray(){
        // return type ## _frompointer((*self).name);
        return type ## _frompointer((*self).name);
    }
    
    %pythoncode %{
       __swig_getmethods__["name"] = name ## _getarray
       if _newclass: name = _swig_property(name ## _getarray, __swig_setmethods__["name"])
    %}
}
%enddef
// /macro definition for referring to array in class like o.a[i] from python



// macro definition for referring to array in class 
// like class itself is array, i. e. o[i] from python
// http://stackoverflow.com/questions/8776328/swig-interfacing-c-library-to-python-creating-iterable-python-data-type-from
// no slicing supported!
// define static int myErr = 0; and %array_class(smth, smthArray); first !
// usage example: ARRAYCLASS(Spam, array, q, float);

%define ARRAYCLASS(cls, name, len, type)

%exception cls::__getitem__ {
  assert(!myErr);
  $action
  if (myErr) {
    myErr = 0; // clear flag for next time
    // You could also check the value in $result, but it''s a PyObject here
    SWIG_exception(SWIG_IndexError, "Can not get cls::name : index out of bounds");
  }
}

%exception cls::__setitem__ {
  assert(!myErr);
  $action
  if (myErr) {
    myErr = 0; // clear flag for next time
    // You could also check the value in $result, but it''s a PyObject here
    SWIG_exception(SWIG_IndexError, "Can not set cls::name : index out of bounds");
  }
}

%extend cls{
    type __getitem__(int index) const {
        if (index >= $self->len) {
            myErr = 1;
            return 0;
        }
        return $self->name[index];
        // return (*self).name[index];
    }
    void __setitem__(int index, type value) const {
        if (index >= $self->len) {
            myErr = 1;
            return;
        }
        (*self).name[index] = value;
    };
    int __len__() {
        return (*self).len;
    }
     
}
%enddef
// /macro definition for referring to array in class like o[i] from python


