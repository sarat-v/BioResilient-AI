import { Routes, Route } from 'react-router-dom'
import { AnimatePresence, motion } from 'framer-motion'
import Layout from '@/components/Layout'
import Dashboard from '@/pages/Dashboard'
import PipelinePage from '@/pages/Pipeline'
import CandidatesPage from '@/pages/Candidates'
import GeneDetail from '@/pages/GeneDetail'
import SpeciesPage from '@/pages/Species'
import ResearchAssistant from '@/pages/ResearchAssistant'

const PAGE_TRANSITION = {
  initial: { opacity: 0, x: 8 },
  animate: { opacity: 1, x: 0 },
  exit: { opacity: 0, x: -8 },
  transition: { duration: 0.2 },
}

function Animated({ children }) {
  return <motion.div {...PAGE_TRANSITION}>{children}</motion.div>
}

export default function App() {
  return (
    <Layout>
      <AnimatePresence mode="wait">
        <Routes>
          <Route path="/" element={<Animated><Dashboard /></Animated>} />
          <Route path="/pipeline" element={<Animated><PipelinePage /></Animated>} />
          <Route path="/candidates" element={<Animated><CandidatesPage /></Animated>} />
          <Route path="/candidates/:id" element={<Animated><GeneDetail /></Animated>} />
          <Route path="/research" element={<Animated><ResearchAssistant /></Animated>} />
          <Route path="/species" element={<Animated><SpeciesPage /></Animated>} />
        </Routes>
      </AnimatePresence>
    </Layout>
  )
}
