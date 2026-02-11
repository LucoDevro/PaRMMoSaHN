import pytest

from package.core import ExampleClass, helper_function


def test_helper_function_returns_expected_message():
    assert helper_function() == "This is a helper function."


@pytest.mark.parametrize("value,expected", [(0, 0), (2, 4), (-3, -6)])
def test_example_class_example_method_doubles_value(value, expected):
    assert ExampleClass(value).example_method() == expected
