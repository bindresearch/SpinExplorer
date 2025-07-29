import pkgutil
import importlib
import inspect
import sys

# Get the current module (this __init__.py)
current_module = sys.modules[__name__]


def auto_import_classes(package):
    for _, module_name, _ in pkgutil.iter_modules(package.__path__):
        full_module_name = f"{package.__name__}.{module_name}"
        module = importlib.import_module(full_module_name)

        for name, obj in inspect.getmembers(module, inspect.isclass):
            if obj.__module__ == full_module_name:
                setattr(current_module, name, obj)


# Auto-import classes from sibling modules
auto_import_classes(sys.modules[__name__])
