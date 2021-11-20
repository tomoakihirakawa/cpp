A2VI = (array)=>{
    let vec = new Module.VectorInteger;
    for(let i=0; i<array.length; i++)
	vec.push_back(array[i]);
    return vec;
};

A2VD = (array)=>{
    let vec = new Module.VectorDouble;
    for(let i=0; i<array.length; i++)
	vec.push_back(array[i]);
    return vec;
};


V2A = (v)=>{
    let s = v.size();
    let array = new Array(s);
    for(let i=0; i<s; i++)
	array[i] = v.get(i);
    return array;
};
/////////////////////////////////////

VVI2A = (v)=>{
    return V2A(Module.FlattenVVI(v));
};

VVD2A = (v)=>{
    console.log(typeof v);
    return V2A(Module.FlattenVVD(v));
};


VVF2A = (v)=>{
    return V2A(Module.FlattenVVF(v));
};

/////////////////////////////////////
V2Float34Array = (v)=>{
    return new Float32Array(V2A(v));
};

V2Unit16Array = (v)=>{
    return new Unit16Array(V2A(v));
};

VVD2Float34Array = (v)=>{
    return new Float32Array(VVD2A(v));
};

VVI2Float16Array = (v)=>{
    return new Float16Array(VVI2A(v));
};
//////////////////////////////////////
setA2V = (vec, array) =>{
    vec.resize(0,0);
    array.forEach(a=>vec.push_back(a));
};

getVec = (vec) =>{
    let ret = new Array(vec.size());
    for(let i=0; i<vec.size(); i++)
	ret[i] = vec.get(i);
    return ret;
};
