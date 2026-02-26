#!/bin/bash
# Update DE-LIMP Docker container with latest code from GitHub
# Usage: bash update_docker.sh

echo "Pulling latest code from GitHub..."
git pull origin main

echo "Rebuilding and restarting container..."
docker compose up -d --build

echo "Done! App running at http://localhost:3838"
echo "Run 'docker compose logs -f' to view logs"
