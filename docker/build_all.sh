#!/usr/bin/env bash
#
# Build all BioResilient Docker images.
#
# Usage:
#   bash docker/build_all.sh                     # build only (local cache)
#   bash docker/build_all.sh --push              # build + push to ECR, no local copy kept
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
echo -e "  Push     : ${BOLD}${PUSH}${RESET}"
echo ""

# Ensure a buildx builder exists (needed for --push without local load)
if $PUSH; then
    docker buildx inspect bioresilient-builder &>/dev/null \
        || docker buildx create --name bioresilient-builder --use
    docker buildx use bioresilient-builder
fi

for img in "${IMAGES[@]}"; do
    banner "Building ${img}"
    local_tag="${REGISTRY}/${img}:${TAG}"
    base_arg="${REGISTRY}/base:${TAG}"

    if $PUSH; then
        docker buildx build \
            --platform linux/amd64 \
            --push \
            --build-arg BASE_IMAGE="${base_arg}" \
            -t "$local_tag" \
            -f "${REPO_ROOT}/docker/${img}/Dockerfile" \
            "${REPO_ROOT}" \
            2>&1 | tail -8
        ok "Built + pushed ${local_tag}"
    else
        docker build \
            --platform linux/amd64 \
            --build-arg BASE_IMAGE="${base_arg}" \
            -t "$local_tag" \
            -f "${REPO_ROOT}/docker/${img}/Dockerfile" \
            "${REPO_ROOT}" \
            2>&1 | tail -5
        ok "Built ${local_tag}"
    fi
done

banner "Done"
echo -e "  Images in ECR:"
for img in "${IMAGES[@]}"; do
    echo -e "    ${REGISTRY}/${img}:${TAG}"
done
echo ""

if ! $PUSH; then
    echo -e "  To push: ${BOLD}REGISTRY=your-ecr-url bash docker/build_all.sh --push${RESET}"
    echo ""
fi
