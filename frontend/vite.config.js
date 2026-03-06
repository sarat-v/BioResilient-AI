import { defineConfig } from 'vite'
import react from '@vitejs/plugin-react'
import path from 'path'

export default defineConfig({
  plugins: [react()],
  resolve: {
    alias: { '@': path.resolve(__dirname, './src') },
  },
  server: {
    proxy: {
      '/candidates': 'http://localhost:8000',
      '/species': 'http://localhost:8000',
      '/scores': 'http://localhost:8000',
      '/pipeline': 'http://localhost:8000',
      '/api': 'http://localhost:8000',
    },
  },
  build: {
    outDir: process.env.VERCEL ? 'dist' : '../api/static/dist',
    emptyOutDir: true,
  },
})
