import { useRelevantArticles } from '../../hooks/useRelevantArticles';
import { Books, LinkSimple, X } from '@phosphor-icons/react';
import { Accordion, Button, ButtonGroup, Row, Col, Badge } from 'react-bootstrap';
import { highlightWords } from '../../utils/utils';
import { useUrlParams } from '../../hooks/useUrlParams';
import useArticleStore from '../../stores/articleStore';
import { useEffect, useMemo } from 'react';
import { capitalizeFirstLetter } from '../../utils/utils';
import { motion } from 'framer-motion';
import { BarLoader } from 'react-spinners';
import ArticleTermInput from './ArticleTermInput';

const RelevantArticles = () => {
  const { relevantArticles, relevantArticlesError, isLoading } = useRelevantArticles();
  const {
    params: { terms, searchBy },
  } = useUrlParams();


  const articleTerms = useArticleStore((state) => state.articleTerms);

  const wordsToHighlight = useMemo(() => {
    return articleTerms.map((item) => item.term);
  }, [relevantArticles]);

  useEffect(() => {
    terms.forEach((term) => {
      useArticleStore.getState().addArticleTerm(capitalizeFirstLetter(term), searchBy, false);
    });
  }, [searchBy, terms]);

  if (relevantArticles.length === 0 && articleTerms.every((term) => !term.isRemovable)) {
    return null;
  }

  return (
    <div className={'mb-5'}>
      <h3 className={'d-flex justify-content-center align-items-center mb-2'}>
        <Books className={'text-secondary'} weight={'light'} />
        <div className={'vr mx-2'} />
        Relevant Articles
      </h3>
      <Row>
        <Col className={'d-flex justify-content-center'}>
          <ArticleTermInput />
        </Col>
      </Row>
      <Row className={'mb-4'}>
        <Col>
          <hr className={'w-25 mx-auto'} />
          <div className={'d-flex justify-content-center flex-wrap'}>
            {articleTerms.map((item) => {
              return (
                <motion.div
                  layout
                  initial={{ opacity: 0 }}
                  animate={{ opacity: 1 }}
                  exit={{ opacity: 0 }}
                  transition={{ type: 'spring', stiffness: 260, damping: 20, delay: 0.1 }}
                  key={item.term + item.type}
                >
                  <ButtonGroup key={item.term + item.type} className={'mx-1 py-1'}>
                    <Button variant={'outline-secondary'} style={{ pointerEvents: 'none' }}>
                      {item.term}
                    </Button>
                    {item.isRemovable && (
                      <Button
                        variant={'outline-secondary'}
                        onClick={() => useArticleStore.getState().removeArticleTerm(item.term, item.type)}
                      >
                        <X />
                      </Button>
                    )}
                  </ButtonGroup>
                </motion.div>
              );
            })}
          </div>
        </Col>
      </Row>
      <Row>
        {isLoading && (
          <div className={'d-flex justify-content-center mb-4'}>
            <BarLoader className={'text-red'} loading={isLoading} speedMultiplier={2}/>
          </div>
        )}
      </Row>
      {(relevantArticles.length > 0 || articleTerms.some((term) => term.isRemovable === true)) && (
        <motion.div layout transition={{ type: 'spring', stiffness: 260, damping: 20, delay: 0.1 }}
        style={{
          filter: isLoading ? 'blur(0.5em) grayscale(1)' : 'none',
          pointerEvents: isLoading ? 'none' : 'auto',
          }}
        >
          {relevantArticles.length > 0 ? (
            <Accordion>
              {relevantArticles.map((article, index) => (
                <Accordion.Item key={index} eventKey={index.toString()}>
                  <Accordion.Header>
                    <div>
                      <div className={'mb-2'}>{article.title}</div>
                      <div className={'d-flex align-items-center small text-secondary'}>
                        {article.authors && article.authors.join(', ')}
                        <div className={'vr mx-2'} />
                        {article.country}
                        <div className={'vr mx-2'} />
                        {article.pm_year}
                      </div>
                    </div>
                  </Accordion.Header>
                  <Accordion.Body>
                    <div>{highlightWords(article.abstract, wordsToHighlight)}</div>
                    <Row className={'mt-3'}>
                      <Col className={'d-flex flex-column'}>
                        <div className={'d-flex align-items-center mb-2'}>
                          <span className={'text-secondary'}>Published in</span>
                          <div className={'vr mx-2'} />
                          {article.venue_title}
                        </div>
                        <div className={'d-flex align-items-center'}>
                          <span className={'text-secondary'}>Publication year</span>
                          <div className={'vr mx-2'} />
                          {article.venue_year}
                        </div>
                        {article.key_words.length !== 0 && (
                          <div className={'d-flex align-items-center flex-wrap mt-2'}>
                            {article.key_words.map((word, index) => (
                              <Badge key={index} bg={'secondary'} className={'me-2 my-1'}>
                                {word}
                              </Badge>
                            ))}
                          </div>
                        )}
                      </Col>
                      <Col className={'d-flex justify-content-end align-items-center'}>
                        <Button href={article.url} target={'_blank'} variant={'outline-primary'}>
                          <div className={'d-flex align-items-center'}>
                            <LinkSimple weight={'light'} />
                            <div className={'vr mx-2'} />
                            Read on PubMed
                          </div>
                        </Button>
                      </Col>
                    </Row>
                  </Accordion.Body>
                </Accordion.Item>
              ))}
            </Accordion>
          ) : (
            <div className={'d-flex justify-content-center'}>
              <div className={'text-secondary'}>No relevant articles found</div>
            </div>
          )}
        </motion.div>
      )}
    </div>
  );
};

export default RelevantArticles;
