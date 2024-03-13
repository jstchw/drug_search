import { Nav, Row, Col } from 'react-bootstrap';
import SearchBox from '../SearchBox/SearchBox';
import Header from '../Header/Header';
import { motion, useInView, AnimatePresence } from 'framer-motion';
import { FC, useState, useEffect, RefObject } from 'react';

interface ScrollbarProps {
  searchBoxSectionRef: RefObject<HTMLDivElement>;
}

const Scrollbar: FC<ScrollbarProps> = ({ searchBoxSectionRef }) => {
  const [showScrollbar, setShowScrollbar] = useState<boolean>(false);
  const isSearchBoxInView = useInView(searchBoxSectionRef);

  useEffect(() => {
    setShowScrollbar(!!!isSearchBoxInView);
  }, [isSearchBoxInView]);

  return (
    <>
    <AnimatePresence>
      {showScrollbar && (
          <motion.div
            initial={{ y: -100 }}
            animate={{ y: 0 }}
            exit={{ y: -100 }}
            transition={{
              type: 'spring',
              stiffness: 260,
              damping: 20,
              duration: 0.8,
              delay: 0.1,
            }}
            className="fixed-top py-4 rounded shadow-lg bg-body"
          >
            <Nav>
              <Row className="w-100 align-items-center">
                <Col className={'ms-4'}>
                  <Header />
                </Col>

                <Col xs={4} className="d-flex justify-content-center">
                  <SearchBox />
                </Col>

                <Col xs={4}></Col>
              </Row>
            </Nav>
          </motion.div>
      )}
    </AnimatePresence>
    </>
  );
};

export default Scrollbar;