test: ## run all tests
	bash run_tests.sh 1

test-fast: ## run tests except test_optimization on all python versions
	bash run_tests.sh 0

bump-major: ## bump major version
	bumpversion major

bump-minor: ## bump minor version
	bumpversion minor

bump-patch: ## bump patch version
	bumpversion patch