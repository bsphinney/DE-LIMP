# Update DE-LIMP Docker container with latest code from GitHub
# Usage: .\update_docker.ps1

Write-Host "Pulling latest code from GitHub..."
git pull origin main

Write-Host "Rebuilding and restarting container..."
docker compose up -d --build

Write-Host "Done! App running at http://localhost:3838"
Write-Host "Run 'docker compose logs -f' to view logs"
