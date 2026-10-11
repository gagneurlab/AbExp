import pickle

import numpy as np
import pytest
from sklearn.exceptions import NotFittedError
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import PolynomialFeatures
from sklearn.utils.validation import check_is_fitted

from abexp.utils.models import AbExpBinaryClassifier

X = np.array([[0.0, 1.0], [1.0, 0.0], [1.0, 1.0], [2.0, 1.0], [0.0, 0.0], [2.0, 2.0]])
Y = np.array([False, False, True, True, False, True])


def test_wrapper_is_not_fitted_before_fit():
    with pytest.raises(NotFittedError):
        check_is_fitted(AbExpBinaryClassifier(LogisticRegression))


def test_wrapper_is_fitted_after_fit():
    check_is_fitted(AbExpBinaryClassifier(LogisticRegression).fit(X, Y))


def test_pipeline_with_the_wrapper_as_last_step_predicts():
    # Pipeline.predict checks that its last step is fitted, and raises NotFittedError since scikit-learn 1.8
    pipeline = Pipeline([
        ("polynomial_features", PolynomialFeatures(interaction_only=True, include_bias=False)),
        ("classifier", AbExpBinaryClassifier(LogisticRegression)),
    ]).fit(X, Y)

    x_poly = PolynomialFeatures(interaction_only=True, include_bias=False).fit_transform(X)
    expected = LogisticRegression().fit(x_poly, Y).predict_proba(x_poly)[:, 1]
    np.testing.assert_array_equal(pipeline.predict(X), expected)


class _PickledByVersion010:
    """Pickles as aberrant-expression 0.1.0 pickled a fitted AbExpBinaryClassifier: its class and its instance dict"""

    def __init__(self, state):
        self.state = state

    def __reduce__(self):
        # pickle.loads creates the AbExpBinaryClassifier without __init__ and then sets the instance dict
        return object.__new__, (AbExpBinaryClassifier,), self.state


def test_wrapper_pickled_by_version_0_1_0_predicts():
    # the predict jobs unpickle models that older releases trained
    fitted_model = LogisticRegression().fit(X, Y)
    pickled = pickle.dumps(_PickledByVersion010({"model_class": LogisticRegression, "model": fitted_model}))

    model = pickle.loads(pickled)

    assert type(model) is AbExpBinaryClassifier
    check_is_fitted(model)
    np.testing.assert_array_equal(model.predict(X), fitted_model.predict_proba(X)[:, 1])
