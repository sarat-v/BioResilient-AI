/** @type {import('tailwindcss').Config} */
export default {
  content: ['./index.html', './src/**/*.{js,ts,jsx,tsx}'],
  theme: {
    extend: {
      colors: {
        // Off-white canvas
        canvas:   '#F7F5F0',
        'canvas-2': '#EDE9E2',
        surface:  '#FFFFFF',
        'surface-2': '#F3F0EC',
        border:   '#E5E0D8',
        'border-2': '#D9D3C9',

        // Text scale (warm, not cool)
        ink:      '#17140F',
        'ink-2':  '#5A5249',
        'ink-3':  '#9E9890',

        // Accent — vivid electric blue
        accent:        '#1B45D4',
        'accent-hover':'#1538B8',
        'accent-light':'#EEF2FF',
        'accent-mid':  '#C7D2FE',

        // Bio green for success/validated
        bio:       '#047857',
        'bio-bg':  '#ECFDF5',
        'bio-ring':'#A7F3D0',

        // Supporting palette
        violet:     '#7C3AED',
        'violet-bg':'#F5F3FF',
        'violet-ring':'#DDD6FE',

        amber:      '#B45309',
        'amber-bg': '#FFFBEB',
        'amber-ring':'#FDE68A',

        slate:      '#475569',
        'slate-bg': '#F8FAFC',
        'slate-ring':'#CBD5E1',

        danger:     '#B91C1C',
        'danger-bg':'#FEF2F2',
        warning:    '#B45309',
        success:    '#047857',
      },
      fontFamily: {
        sans: ['"Plus Jakarta Sans"', 'system-ui', 'sans-serif'],
        mono: ['"JetBrains Mono"', 'monospace'],
      },
      backgroundImage: {
        'gradient-accent': 'linear-gradient(135deg, #1B45D4, #7C3AED)',
        'gradient-bio':    'linear-gradient(135deg, #059669, #1B45D4)',
        'gradient-warm':   'linear-gradient(135deg, #F7F5F0, #EDE9E2)',
      },
      animation: {
        'fade-in':   'fadeIn 0.4s ease-out',
        'slide-up':  'slideUp 0.4s ease-out',
        'pulse-slow':'pulse 3s cubic-bezier(0.4, 0, 0.6, 1) infinite',
      },
      keyframes: {
        fadeIn:  { from: { opacity: '0' }, to: { opacity: '1' } },
        slideUp: { from: { opacity: '0', transform: 'translateY(12px)' }, to: { opacity: '1', transform: 'translateY(0)' } },
      },
      boxShadow: {
        'card':    '0 1px 3px rgba(23, 20, 15, 0.06), 0 1px 2px rgba(23, 20, 15, 0.04)',
        'card-md': '0 4px 16px rgba(23, 20, 15, 0.08), 0 1px 4px rgba(23, 20, 15, 0.04)',
        'card-lg': '0 8px 32px rgba(23, 20, 15, 0.10), 0 2px 8px rgba(23, 20, 15, 0.06)',
        'accent':  '0 4px 14px rgba(27, 69, 212, 0.25)',
        'focus':   '0 0 0 3px rgba(27, 69, 212, 0.20)',
      },
      borderRadius: {
        '2xl': '1rem',
        '3xl': '1.5rem',
      },
    },
  },
  plugins: [],
}
