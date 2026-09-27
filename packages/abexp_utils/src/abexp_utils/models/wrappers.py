import sys

import numpy as np
from sklearn.base import BaseEstimator
from sklearn.base import RegressorMixin

__all__ = [
    "ModelWrapper",
    "AbExpBinaryClassifier",
    "AbExpZscoreRegressor",
    "NullRegressor",
]


def _import_string(dotted_path):
    """
    Import a dotted module path and return the attribute/class designated by the
    last name in the path. Raise ImportError if the import failed.
    """
    try:
        module_path, class_name = dotted_path.rsplit('.', 1)
    except ValueError:
        msg = "%s doesn't look like a module path" % dotted_path
        import six
        six.reraise(ImportError, ImportError(msg), sys.exc_info()[2])

    from importlib import import_module
    module = import_module(module_path)

    try:
        return getattr(module, class_name)
    except AttributeError:
        msg = 'Module "%s" does not define a "%s" attribute/class' % (
            module_path, class_name)
        import six
        six.reraise(ImportError, ImportError(msg), sys.exc_info()[2])


class ModelWrapper(BaseEstimator):

    def __init__(self, model_class, model_instance=None, **kwargs):
        # convert string like "sklearn.linear_model.LinearRegression" to actual class
        if isinstance(model_class, str):
            model_class = _import_string(model_class)

        self.model_class = model_class
        if model_instance is None:
            self.model = model_class(**kwargs)
        else:
            self.model = model_instance

    def fit(self, *args, **kwargs):
        fitted_model = self.model.fit(*args, **kwargs)
        if fitted_model is not None and fitted_model is not self.model:
            self.model = fitted_model
        return self

    def get_params(self, *args, **kwargs):
        retval = self.model.get_params()
        retval["model_class"] = self.model_class
        return retval

    def __getattr__(self, *args):
        return self.model.__getattribute__(*args)

    def __repr__(self, *args, **kwargs):
        return f"{self.__class__.__name__}({repr(self.model, *args, **kwargs)})"

    def __str__(self, *args, **kwargs):
        return f"{self.__class__.__name__}({str(self.model, *args, **kwargs)})"


class AbExpBinaryClassifier(ModelWrapper):

    def predict(self, *args, **kwargs):
        return self.predict_proba(*args, **kwargs)

    def predict_proba(self, *args, **kwargs):
        return self.model.predict_proba(*args, **kwargs)[:, np.argwhere(self.model.classes_).item()]


class AbExpZscoreRegressor(ModelWrapper):

    def predict_proba(self, *args, **kwargs):
        return -self.model.predict(*args, **kwargs)


class NullRegressor(RegressorMixin):
    """
    Dummy model that returns the mean of the features as output.
    """
    def __init__(self):
        pass

    def fit(self, X=None, y=None):
        return self

    def predict(self, X=None):
        X = np.asanyarray(X)
        return np.mean(X, axis=-1)
    
    def predict_proba(self, *args, **kwargs):
        return self.predict(*args, **kwargs)

    def get_params(self, deep = False):
        return {}


