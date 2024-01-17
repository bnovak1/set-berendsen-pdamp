dist: ## builds source and wheel package 
	python setup.py sdist bdist_wheel 

clean-build: ## remove build artifacts
	rm -fr build/
	rm -fr dist/
	rm -fr .eggs/

test: ## run all tests on all python & LAMMPS versions
	bash run_tests.sh 1

test-fast: ## run tests except test_optimization on all python versions
	bash run_tests.sh 0

release: dist ## package and upload a release
	twine upload dist/*

bump-major: ## bump major version
	bumpversion major

bump-minor: ## bump minor version
	bumpversion minor

bump-patch: ## bump patch version
	bumpversion patch