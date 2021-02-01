
# Variables

makefile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
makefile_dir := $(dir $(makefile_path))

venv_path := ${makefile_dir}.venv
env_path = $(makefile_dir)docker/.env
venv_ok_path := ${venv_path}/venv_ok.txt
docker_compose_file = ${makefile_dir}docker/docker-compose.yml

# helpers for notebooks

examples_dir = ${makefile_dir}examples
notebooks_dir = ${makefile_dir}notebooks
project_dir = ${makefile_dir}

examples_dir = ${makefile_dir}examples
notebooks_config_path = ${notebooks_dir}/config.ini
project_config_path = ${makefile_dir}config.ini
config_txt = project_dir=${project_dir}

ifdef OS
	VENV_CREATE_PYTHON = py -3
	DEL_CMD = del
	THEN_RUN = &
	pip_bin = ${venv_path}/Scripts/pip.exe
	precommit_bin = ${venv_path}/Scripts/pre-commit.exe
else
	VENV_CREATE_PYTHON = python3
	DEL_CMD = rm -f
	THEN_RUN = ;
	pip_bin = ${venv_path}/bin/pip
	precommit_bin = ${venv_path}/bin/pre-commit
endif

# Build

${env_path}:
	$(info Writing ${env_path})
	$(file > ${env_path}, DISPLAY=${DISPLAY})

env: ${env_path}

${venv_ok_path}:
	${VENV_CREATE_PYTHON} -m venv --clear "${venv_path}"
	"${pip_bin}" install pre-commit==2.8.2
	"${precommit_bin}" install
	echo ok > "${venv_ok_path}"

venv: ${venv_ok_path}

build: env venv
	docker-compose -f "${docker_compose_file}" build --no-cache

# Run

.PHONY: run-fenics, up

run-fenics:
	docker-compose -f "${docker_compose_file}" run fenics

up:
	docker-compose -f "${docker_compose_file}" up

.PHONY: tests, tests-docker

tests:
	pytest ./tests ./pkgs/fenics_utils/tests ./pkgs/paraview_utils/tests
	
tests-docker:
	docker-compose -f "${docker_compose_file}" run fenics "pytest ./tests ./pkgs/fenics_utils/tests ./pkgs/paraview_utils/tests"


# Notebook config

${examples_config_path}:
	$(file > ${examples_config_path}, [Config])
	$(file >> ${examples_config_path}, ${config_txt})

${notebooks_config_path}:
	$(file > ${notebooks_config_path}, [Config])
	$(file >> ${notebooks_config_path}, ${config_txt})

${project_config_path}:
	$(file > ${project_config_path}, [Config])
	$(file >> ${project_config_path}, ${config_txt})

examples_config: ${examples_config_path}
notebooks_config: ${notebooks_config_path}
project_config: ${project_config_path}

config: examples_config notebooks_config project_config

notebooks: config
	jupytext --to notebook ${examples_dir}/twoD_method_of_images.py -o ${notebooks_dir}/twoD_method_of_images.ipynb

server: config notebooks
	jupyter lab --ip=0.0.0.0 --port=8888 --allow-root
