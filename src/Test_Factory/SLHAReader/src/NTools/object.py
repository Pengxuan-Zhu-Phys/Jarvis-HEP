#!/usr/bin/env python3

STYLE = {
        'fore':
        {   # 前景色
            'black'    : 30,   #  黑色
            'red'      : 31,   #  红色
            'green'    : 32,   #  绿色
            'yellow'   : 33,   #  黄色
            'blue'     : 34,   #  蓝色
            'purple'   : 35,   #  紫红色
            'cyan'     : 36,   #  青蓝色
            'white'    : 37,   #  白色
        },

        'back' :
        {   # 背景
            'black'     : 40,  #  黑色
            'red'       : 41,  #  红色
            'green'     : 42,  #  绿色
            'yellow'    : 43,  #  黄色
            'blue'      : 44,  #  蓝色
            'purple'    : 45,  #  紫红色
            'cyan'      : 46,  #  青蓝色
            'white'     : 47,  #  白色
        },

        'mode' :
        {   # 显示模式
            'mormal'    : 0,   #  终端默认设置
            'bold'      : 1,   #  高亮显示
            'underline' : 4,   #  使用下划线
            'blink'     : 5,   #  闪烁
            'invert'    : 7,   #  反白显示
            'hide'      : 8,   #  不可见
        },

        'default' :
        {
            'end' : 0,
        },
}

def AcquireAttr(obj,attr,value):
    if not hasattr(obj,attr):
        setattr(obj,attr,value)
    return getattr(obj,attr)

class lazyproperty: # P271 in Python Cookbook
    def __init__(self,func):
        self.func=func
    def __get__(self,instance,cls):
        if instance is None:
            return self
        else:
            value=self.func(instance)
            setattr(instance,self.func.__name__,value)
            return value

class lazyclass: # P271 in Python Cookbook
    def __init__(self,func,arg):
        self.func=func
        self.args=args
    def __get__(self,instance,cls):
        if instance is None:
            return self
        else:
            value=self.func(instance,*self.args)
            setattr(instance,self.func.__name__,value)
            return value

def MergeAttributes(*objects):
    new_obj=capsule()
    for obj in objects:
        for key,value in obj.__dict__.items():
            setattr(new_obj,key,value)
    return new_obj

def Error(*arguments,**keywords):
    ColorPrint(1,31,'',"\nError:\n  ",end='')
    ColorPrint(0,31,'',*arguments,**keywords)
    exit()


class capsule():
    '''container to store data in it,
    thus data can be accessed by .__dict__
    '''
    def Merge(self,*others):
        self.__dict__= MergeAttributes(self,*others).__dict__
    def MoveTo(self,destinations):
        keys=destinations.keys()
        if self.documents.keys()!=keys:
            Error('keys do not match while moving files')
        else:
            for key in keys:
                shutil.move(self.documents[key],destinations[key])
            self.documents.update(destinations)
    def CopyTo(self,destinations):
        keys=destinations.keys()
        for key in keys:
            try:
                shutil.copy(self.documents[key],destinations[key])
            except TypeError:
                if self.documents[key] is None:
                    destinations[key]=None
                else:
                    raise                
        self.documents.update(destinations)
    def EqualTo(self,other,block_list=None):
        if block_list is None:
            block_list=list(self.__dict__.keys())
            #print(block_list)
        eq=all([
            getattr(self,block)==getattr(other,block,list())
            for block in block_list
        ])
        return eq
    def __eq__(self,other):
        return self.EqualTo(other)


def ColorPrint(mode='',fore='',back='',*arguments,**keywords):
    if mode in STYLE['mode'].keys():
        mode='%s' %STYLE['mode'][mode]
    elif mode in STYLE['mode'].values():
        mode='%s' %mode
    else:
        mode=''
    
    if fore in STYLE['fore'].keys():
        fore='%s' %STYLE['fore'][fore]
    elif fore in STYLE['fore'].values():
        fore='%s' %fore
    else:
        fore=''

    if back in STYLE['back'].keys():
        back='%s' %STYLE['back'][back]
    elif back in STYLE['back'].values():
        back='%s' %back
    else:
        back=''
    

    # mode  = '%s' % STYLE['mode'][mode] if mode in STYLE['mode'].keys() else ''
    # fore  = '%s' % STYLE['fore'][fore] if fore in STYLE['fore'].keys() else ''
    # back  = '%s' % STYLE['back'][back] if back in STYLE['back'].keys() else ''    

    style = ';'.join([s for s in [mode, fore, back] if s])
    begin = '\033[%sm' % style if style else ''
    end   = '\033[%sm' % STYLE['default']['end'] if style else ''
    print(begin,end='')
    print(*arguments,**keywords)
    print(end,end='')