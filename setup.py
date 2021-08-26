
from setuptools import setup


setup(
    name="Plasmaman",
    version=open("VERSION.txt").read().strip(),
    author="Evgeny P. Kurbatov",
    author_email="kurbatov@inasan.ru",
    packages=['plasmaman'],
    url='https://github.com/evgenykurbatov/plasmaman.git',
    description="Primitives for modeling the interaction of plasma and electromagnetic fields",
    long_description=open("README.md").read(),
    install_requires=[
        "NumPy",
        "SciPy"
    ],
)
