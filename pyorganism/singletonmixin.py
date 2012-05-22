"""
A Python Singleton mixin class that makes use of some of the ideas
found at http://c2.com/cgi/wiki?PythonSingleton. Just inherit
from it and you have a singleton. No code is required in
subclasses to create singleton behavior -- inheritance from
Singleton is all that is needed.

Singleton creation is threadsafe.

USAGE:

Just inherit from Singleton. If you need a constructor, include
an __init__() method in your class as you usually would. However,
if your class is S, you instantiate the singleton using S.get_instance()
instead of S(). Repeated calls to S.get_instance() return the
originally-created instance.

For example:


class S(Singleton):

    def __init__(self, a, b=1):
        pass

S1 = S.get_instance(1, b=3)


Most of the time, that"s all you need to know. However, there are some
other useful behaviors. Read on for a full description:

1) Getting the singleton:

    S.get_instance()

returns the instance of S. If none exists, it is created.

3) Use __S.__init__() for instantiation processing,
since S.get_instance() runs S.__init__(), passing it the args it has received.

If no data needs to be passed in at instantiation time, you don"t need S.__init__().

4) If S.__init__(.) requires parameters, include them ONLY in the
first call to S.get_instance(). If subsequent calls have arguments,
a SingletonException is raised by default.

If you find it more convenient for subsequent calls to be allowed to
have arguments, but for those argumentsto be ignored, just include
"ignore_subsequent = True" in your class definition, i.e.:

  class S(Singleton):

      ignore_subsequent = True

      def __init__(self, a, b=1):
          pass

5) For testing, it is sometimes convenient for all existing singleton
instances to be forgotten, so that new instantiations can occur. For that
reason, a forget_all_singletons() function is included. Just call

  forget_all_singletons()

and it is as if no earlier instantiations have occurred.

6) As an implementation detail, classes that inherit
from Singleton may not have their own __new__
methods. To make sure this requirement is followed,
an exception is raised if a Singleton subclass includ
es __new__. This happens at subclass instantiation
time (by means of the MetaSingleton metaclass.


By Gary Robinson, grobinson@flyfi.com. No rights reserved --
placed in the public domain -- which is only reasonable considering
how much it owes to other people"s code and ideas which are in the
public domain. The idea of using a metaclass came from
a  comment on Gary"s blog (see
http://www.garyrobinson.net/2004/03/python_singleto.html#comments).
Other improvements came from comments and email from other
people who saw it online. (See the blog post and comments
for further credits.)

Not guaranteed to be fit for any particular purpose. Use at your
own risk.
"""

import threading

class SingletonException(StandardError):
    pass

_st_singletons = set()
_lock_for_singletons = threading.RLock()
_lock_for_singleton_creation = threading.RLock()   # Ensure only one instance of each Singleton
                                                # class is created.  This is not bound to the 
                                                # individual Singleton class since we need to
                                                # ensure that there is only one mutex for each
                                                # Singleton class, which would require having
                                                # a lock when setting up the Singleton class,
                                                # which is what this is anyway.  So, when any
                                                # Singleton is created, we lock this lock and
                                                # then we don"t need to lock it again for that
                                                # class.

def _create_singleton_instance(cls, args, kw_args):
    _lock_for_singleton_creation.acquire()
    try:
        if cls._is_instantiated(): # some other thread got here first
            return

        instance = cls.__new__(cls)
        try:
            instance.__init__(*args, **kw_args)
        except TypeError as e:
            if e.message.find("__init__() takes") != -1:
                raise SingletonException("If the singleton requires __init__ "\
                        "args, supply them on only on the first call to "\
                        "get_instance().")
            else:
                raise
        cls.c_instance = instance
        _add_singleton(cls)
    finally:
        _lock_for_singleton_creation.release()

def _add_singleton(cls):
    _lock_for_singletons.acquire()
    try:
        assert cls not in _st_singletons
        _st_singletons.add(cls)
    finally:
        _lock_for_singletons.release()

def _remove_singleton(cls):
    _lock_for_singletons.acquire()
    try:
        if cls in _st_singletons:
            _st_singletons.remove(cls)
    finally:
        _lock_for_singletons.release()

def forget_all_singletons():
    """This is useful in tests, since it is hard to know which singletons need
    to be cleared to make a test work."""
    _lock_for_singletons.acquire()
    try:
        for cls in _st_singletons.copy():
            cls._forget_class_instance_reference_for_testing()

        # Might have created some Singletons in the process of tearing down.
        # Try one more time - there should be a limit to this.
        i_num_singletons = len(_st_singletons)
        if len(_st_singletons) > 0:
            for cls in _st_singletons.copy():
                cls._forget_class_instance_reference_for_testing()
                i_num_singletons -= 1
                assert i_num_singletons == len(_st_singletons),\
                        "Added a singleton while destroying " + str(cls)
        assert len(_st_singletons) == 0, _st_singletons
    finally:
        _lock_for_singletons.release()

class MetaSingleton(type):
    def __new__(metaclass, name, bases, dct):
        if dct.has_key("__new__"):
            raise SingletonException("Can not override __new__ in a Singleton")
        return super(MetaSingleton, metaclass).__new__(metaclass, name, bases, dct)

    def __call__(cls, *args, **kw_args):
        raise SingletonException("Singletons may only be instantiated through"\
                " get_instance()")

class Singleton(object):
    __metaclass__ = MetaSingleton

    def get_instance(cls, *args, **kw_args):
        """
        Call this to instantiate an instance or retrieve the existing instance.
        If the singleton requires args to be instantiated, include them the first
        time you call get_instance.
        """
        if cls._is_instantiated():
            if (args or kw_args) and not hasattr(cls, "ignore_subsequent"):
                raise SingletonException("Singleton already instantiated, but"\
                        " get_instance() called with args.")
        else:
            _create_singleton_instance(cls, args, kw_args)
        return cls.c_instance
    get_instance = classmethod(get_instance)

    def _is_instantiated(cls):
        # Don"t use hasattr(cls, "c_instance"), because that screws things up if there is a singleton that
        # extends another singleton.  hasattr looks in the base class if it doesn"t find in subclass.
        return "c_instance" in cls.__dict__
    _is_instantiated = classmethod(_is_instantiated)

    # This can be handy for public use also
    is_instantiated = _is_instantiated

    def _forget_class_instance_reference_for_testing(cls):
        """
        This is designed for convenience in testing -- sometimes you
        want to get rid of a singleton during test code to see what
        happens when you call get_instance() under a new situation.

        To really delete the object, all external references to it
        also need to be deleted.
        """
        try:
            if hasattr(cls.c_instance, "_prepare_to_forget_singleton"):
                # tell instance to release anything it might be holding onto.
                cls.c_instance._prepare_to_forget_singleton()
            del cls.c_instance
            _remove_singleton(cls)
        except AttributeError:
            # run up the chain of base classes until we find the one that has the instance
            # and then delete it there
            for base_cls in cls.__bases__: 
                if issubclass(base_cls, Singleton):
                    base_cls._forget_class_instance_reference_for_testing()
    _forget_class_instance_reference_for_testing = classmethod(_forget_class_instance_reference_for_testing)

