#!/usr/bin/env bash
#
# Build all BioResilient Docker images.
#
# Usage:
#   bash docker/build_all.sh                     # build only
#   bash docker/build_all.sh --push              # build + push to registry
#   REGISTRY=xxxx.dkr.ecr.ap-south-1.amazonaws.com bash docker/build_all.sh --push
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REGISTRY="${REGISTRY:-bioresilient}"
TAG="${TAG:-latest}"
PUSH=false
[[ "${1:-}" == "--push" ]] && PUSH=true

IMAGES=(base orthofinder align phylo hyphy esm clinical)

GREEN="\033[0;32m"; CYAN="\033[0;36m"; BOLD="\033[1m"; RESET="\033[0m"
banner() { echo -e "\n${BOLD}${CYAN}━━━  $1  ━━━${RESET}\n"; }
ok()     { echo -e "${GREEN}  ✓ $1${RESET}"; }

banner "Building BioResilient Docker Images"
echo -e "  Registry : ${BOLD}${REGISTRY}${RESET}"
echo -e "  Tag      : ${BOLD}${TAG}${RESET}"
echo ""

for img in "${IMAGES[@]}"; do
    banner "Building ${img}"
    local_tag="${REGISTRY}/${img}:${TAG}"

    docker build \
        -t "$local_tag" \
        -f "${REPO_ROOT}/docker/${img}/Dockerfile" \
        "${REPO_ROOT}" \
        2>&1 | tail -5

    ok "Built ${local_tag}"

    if $PUSH; then
        echo "  Pushing ${local_tag} ..."
        docker push "$local_tag"
        ok "Pushed ${local_tag}"
    fi
done

banner "All images built"
echo -e "  Images:"
for img in "${IMAGES[@]}"; do
    echo -e "    ${REGISTRY}/${img}:${TAG}"
done
echo ""

if ! $PUSH; then
    echo -e "  To push: ${BOLD}REGISTRY=your-ecr-url bash docker/build_all.sh --push${RESET}"
    echo ""
fi
