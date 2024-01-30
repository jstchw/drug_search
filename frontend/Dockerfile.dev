FROM node:21-slim

WORKDIR /usr/src/app

COPY package.json yarn.lock ./

RUN yarn install

COPY src ./src
COPY public ./public
COPY tsconfig.json ./

EXPOSE 3000

ENV CHOKIDAR_USEPOLLING=true

CMD ["yarn", "start"]