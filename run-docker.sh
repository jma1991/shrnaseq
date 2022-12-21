export TEST_ENV=false; \
docker compose -f docker/docker-compose.yml --env-file .env down -v; \
docker compose -f docker/docker-compose.yml --env-file .env up --build; \