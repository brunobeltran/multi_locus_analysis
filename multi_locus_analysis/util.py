import warnings

code_only_warning = "This function is only kept around 'just in case' I end "
"up needing it. While it was likely tested at some point, it was never "
"intended for general use and likely only works in the very narrow scope "
"for which it was written. If you want to accomplish a similar goal, you "
"will likely have to read through the code itself and make significant "
"changes. To ignore this warning, pass `ignore_bruno=True`."


def mark_code_only(f):
    def f_wrap(*args, **kwargs):
        if kwargs.pop('ignore_bruno', None) is not None:
            warnings.warn(code_only_warning)
