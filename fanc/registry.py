class_name_dict = {}
class_id_dict = {}


def get_class_by_name(classname):
    """
    Get the node class matching the `classname`.
    If the name is not registered, a ``TypeError`` is raised.
    """

    # Get the class object corresponding to `classname`.
    if classname not in class_name_dict:
        raise TypeError("there is no registered node class named ``{}``".format(classname))

    return class_name_dict[classname]
